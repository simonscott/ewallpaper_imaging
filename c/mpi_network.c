#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <pthread.h>
#include <string.h>
#include <math.h>
#include "virtual_cpu.h"
#include "network.h"

// MPI Network
int proc_x, proc_y;
int nx_proc, ny_proc;

// Network Parameters
int num_x, num_y;
int nx_cpu, ny_cpu;
cpu* cpus;

//================================================================================
//============================ Utilities =========================================
//================================================================================

int is_local(int x, int y){
  return
    x >= proc_x * nx_cpu && x < (proc_x+1)*nx_cpu &&
    y >= proc_y * ny_cpu && y < (proc_y+1)*ny_cpu;
}

void to_local_coords(int* x, int* y){
  if(!is_local(*x, *y)){
    printf("CPU %d,%d is not local to node %d,%d.\n", *x, *y, proc_x, proc_y);
    exit(-1);
  }
  
  x[0] %= nx_cpu;
  y[0] %= ny_cpu;
}

int get_rank(int proc_x, int proc_y){
  return proc_x + proc_y * nx_proc;
}

int pack_tag(int x, int y, int dir){
  return (dir << 16) + (x << 8) + y;
}

void unpack_tag(int tag, int* x, int* y, int* dir){
  *dir = tag >> 16;
  *x = (tag >> 8) & 0xFF;
  *y = tag & 0xFF;
}

cpu* get_cpu(int x, int y){
  to_local_coords(&x, &y);
  return &cpus[x + y * nx_cpu];
}

//================================================================================
//========================== Network Lock ========================================
//================================================================================

pthread_cond_t network_signal;
pthread_mutex_t network_lock;

void notify_network_listeners(int* status_bit){
  pthread_mutex_lock(&network_lock);
  status_bit[0] = 1;
  pthread_cond_broadcast(&network_signal);
  pthread_mutex_unlock(&network_lock);
}

void wait_for_network(int* status_bit){
  pthread_mutex_lock(&network_lock);
  while(!status_bit[0])
     pthread_cond_wait(&network_signal, &network_lock);
  pthread_mutex_unlock(&network_lock);
}

//================================================================================
//=========================== MPI Thread =========================================
//================================================================================

//--------------------------------------------------------------------------------
//--------------------------- Magic Constants ------------------------------------
const int ack_tag = 0x4FFFF;

//--------------------------------------------------------------------------------
//--------------------------- Status Bits ----------------------------------------
int* send_status;

void init_status_bits(){
  int max_cpu = nx_cpu > ny_cpu ? nx_cpu : ny_cpu;
  send_status = (int*)malloc(4 * max_cpu * sizeof(int));
}

int* get_status_bit(int dir, int x, int y){
  // Check if x,y is local
  if(!is_local(x,y)){
    printf("Cannot retrieve status bit for non-local cpu %d, %d.\n", x, y);
    exit(-1);
  }
  
  to_local_coords(&x,&y);
  int max_cpu = nx_cpu > ny_cpu ? nx_cpu : ny_cpu;
  if(dir == up || dir == down)
    return &(send_status[x + dir * max_cpu]);
  else if(dir == left || dir == right)
    return &(send_status[y + dir * max_cpu]);
  else{
    printf("Unrecognized direction: %d\n", dir);
    exit(-1);
  }
}

//--------------------------------------------------------------------------------
//------------------------------ Message Queue -----------------------------------
typedef struct {
  int x;
  int y;
  int dir;
  int size;
} message;

typedef struct {
  char* buffer;
  char* buffer_top;
  message* items;
  int size;
  int capacity;
} queue;

void init_queue(queue* q, int capacity, int message_size){
  q->buffer = (char*)malloc(capacity * message_size);
  q->buffer_top = q->buffer;
  q->items = (message*)malloc(capacity * sizeof(message));
  q->size = 0;
  q->capacity = capacity;
}

void add_to_queue(queue* q, message msg, void* msg_buffer){
  //Check for space
  if(q->size == q->capacity){
    printf("Queue is full.\n");
    exit(-1);
  }

  //Copy message into buffer
  if(msg_buffer)
    memcpy(q->buffer_top, msg_buffer, msg.size);
  q->buffer_top += msg.size;
  //Store message
  q->items[q->size] = msg;
  (q->size)++;
}

char* get_message_buffer(queue* q, int i){
  //Check Bounds
  if(i < 0 || i >= q->size){
    printf("Queue index out of bounds %d\n", i);
    exit(-1);
  }

  char* buffer_start = q->buffer;
  for(int k=0; k<i; k++)
    buffer_start += q->items[k].size;
  return buffer_start;
}

void remove_from_queue(queue* q, int i){
  //Check Bounds
  if(i < 0 || i >= q->size){
    printf("Queue index out of bounds %d\n", i);
    exit(-1);
  }

  // Get message i
  message msg = q->items[i];
  // Compute msg buffer index
  char* buffer_start = get_message_buffer(q, i);
  // Move messages over
  memmove(buffer_start, buffer_start + msg.size,
          q->buffer_top - (buffer_start + msg.size));
  q->buffer_top -= msg.size;
  // Shift messages down
  for(int k=i; k<q->size-1; k++)
    q->items[k] = q->items[k+1];
  (q->size)--;  
}

message peek_queue(queue* q){
  //Check non empty
  if(q->size == 0){
    printf("Queue is empty.\n");
    exit(-1);
  }
  return q->items[0];
}

void pop_queue(queue* q){
  remove_from_queue(q, 0);
}

//--------------------------------------------------------------------------------
//---------------------------- Send and Receive Queues ---------------------------

queue* send_queue;
queue* receive_queue;

void init_network_queues(){
  //Network size
  int queue_size = 16 * 8;
  int msg_size = 256 * 2 * sizeof(int);
  
  send_queue = (queue*)malloc(sizeof(queue));
  receive_queue = (queue*)malloc(sizeof(queue));
  init_queue(send_queue, queue_size, msg_size);
  init_queue(receive_queue, queue_size, msg_size);
}

void add_to_receive_queue(message msg, void* msg_buffer){
  pthread_mutex_lock(&network_lock);
  add_to_queue(receive_queue, msg, msg_buffer);
  pthread_mutex_unlock(&network_lock);
}

void add_to_send_queue(message msg, void* msg_buffer){
  pthread_mutex_lock(&network_lock);
  add_to_queue(send_queue, msg, msg_buffer);
  pthread_mutex_unlock(&network_lock);
}

void pop_send_queue(){
  pthread_mutex_lock(&network_lock);
  pop_queue(send_queue);
  pthread_mutex_unlock(&network_lock);
}

void pop_receive_queue(){
  pthread_mutex_lock(&network_lock);
  pop_queue(receive_queue);
  pthread_mutex_unlock(&network_lock);
}

void remove_from_receive_queue(int i){
  pthread_mutex_lock(&network_lock);
  remove_from_queue(receive_queue, i);
  pthread_mutex_unlock(&network_lock);  
}

void process_receive_queue(){
  int index = 0;
  while(index < receive_queue->size){
    message msg = receive_queue->items[index];
    
    // Compute receiving port
    cpu* receiving_cpu;
    if(msg.dir == up)
      receiving_cpu = get_cpu(msg.x, proc_y * ny_cpu);
    else if(msg.dir == down)
      receiving_cpu = get_cpu(msg.x, (proc_y + 1) * ny_cpu - 1);
    else if(msg.dir == left)
      receiving_cpu = get_cpu((proc_x + 1) * nx_cpu - 1, msg.y);
    else if(msg.dir == right)
      receiving_cpu = get_cpu(proc_x * nx_cpu, msg.y);
    cpu_port* port = get_port(receiving_cpu, opposite_dir(msg.dir));

    // If port is not available then go to the next message
    if(!port->space_available){
      index++;
    }
    // Otherwise process the message
    else{
      // Send port message, and remove from queue
      send_port(port, get_message_buffer(receive_queue, index), msg.size);
      remove_from_receive_queue(index);

      // Send acknowledgement message
      msg.size = 0;
      add_to_send_queue(msg, NULL);
    }
  }
}

//--------------------------------------------------------------------------------
//-------------------------- MPI Thread ------------------------------------------
pthread_t mpi_receive_thread;
volatile int mpi_thread_running;
int max_size;
char* mpi_receive_buffer;
char* mpi_send_buffer;

void init_mpi_thread(){
  // Initialize send_status
  init_status_bits();

  // Initialize network signal and lock
  pthread_cond_init(&network_signal, NULL);
  pthread_mutex_init(&network_lock, NULL);

  // Initialize queues
  init_network_queues();

  // Set to running
  mpi_thread_running = 1;

  // Initialize Buffers
  max_size = 1024 * sizeof(int);
  mpi_receive_buffer = (char*)malloc(max_size);
  mpi_send_buffer = (char*)malloc(max_size);
}

void* mpi_thread(void* args){
  // Initialize Buffers and Status
  MPI_Request receive_request = NULL;
  MPI_Request send_request = NULL;
  MPI_Status receive_status;
  MPI_Status send_status;

  // Initial Receive
  MPI_Irecv(mpi_receive_buffer, max_size, MPI_BYTE,
            MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &receive_request);

  while(mpi_thread_running || send_queue->size > 0){    
    // Check if receive complete
    int receive_buffer_free;
    MPI_Test(&receive_request, &receive_buffer_free, &receive_status);
    if(receive_buffer_free){
      // Acknowledgement message
      if(receive_status.MPI_TAG == ack_tag) {
        message msg = *(message*)mpi_receive_buffer;
        int* status_bit = get_status_bit(msg.dir, msg.x, msg.y);
        notify_network_listeners(status_bit);                
      }
      // Normal message
      else{
        message msg;
        MPI_Get_count(&receive_status, MPI_BYTE, &(msg.size));
        unpack_tag(receive_status.MPI_TAG, &(msg.x), &(msg.y), &(msg.dir));
        add_to_receive_queue(msg, mpi_receive_buffer);
      }
      // Receive again
      MPI_Irecv(mpi_receive_buffer, max_size, MPI_BYTE,
                MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &receive_request);
    }
    
    // Check if there's anything to send    
    if(send_queue->size > 0){
      // Check if the send complete
      int send_buffer_free;
      if(send_request)
        MPI_Test(&send_request, &send_buffer_free, &send_status);
      else
        send_buffer_free = 1;

      // If the buffer is free, then send the message
      if(send_buffer_free){
        // Get first message in queue
        message msg = peek_queue(send_queue);

        // If Acknowledgement message (size == 0)
        if(msg.size == 0){
          // Compute destination rank
          int dest_proc_x = proc_x;
          int dest_proc_y = proc_y;
          step_in_dir(opposite_dir(msg.dir), &dest_proc_x, &dest_proc_y);
          int dest_rank = get_rank(dest_proc_x, dest_proc_y);

          memcpy(mpi_send_buffer, send_queue->items, sizeof(message));
          MPI_Isend(mpi_send_buffer, sizeof(message), MPI_BYTE, dest_rank, ack_tag,
                    MPI_COMM_WORLD, &send_request);
        }
        // Normal Message
        else {
          // Compute destination rank
          int dest_proc_x = proc_x;
          int dest_proc_y = proc_y;
          step_in_dir(msg.dir, &dest_proc_x, &dest_proc_y);
          int dest_rank = get_rank(dest_proc_x, dest_proc_y);
        
          // Compute message tag
          int dest_tag = pack_tag(msg.x, msg.y, msg.dir);

          // Copy message into buffer and send
          memcpy(mpi_send_buffer, send_queue->buffer, msg.size);
          MPI_Isend(mpi_send_buffer, msg.size, MPI_BYTE, dest_rank, dest_tag,
                    MPI_COMM_WORLD, &send_request);
        }
        
        // Pop the message from queue now that it is sent.
        pop_send_queue();
      }
    }

    // Process received messages
    process_receive_queue();
  }
}

void start_mpi_thread(){
  pthread_create(&mpi_receive_thread, NULL, &mpi_thread, NULL);
}

void end_mpi_thread(){
  mpi_thread_running = 0;
  pthread_join(mpi_receive_thread, NULL);
}

//================================================================================
//=================== Send and Receive ===========================================
//================================================================================

void send_message(int x, int y, int direction, void* msg_buffer, int message_size){
  // Compute local coords
  int dest_x = x;
  int dest_y = y;
  step_in_dir(direction, &dest_x, &dest_y);
  // Send to local processor
  if(is_local(dest_x, dest_y)){    
    cpu* dest_proc = get_cpu(dest_x, dest_y);
    cpu_port* port = get_port(dest_proc, opposite_dir(direction));
    send_port(port, msg_buffer, message_size);
  }
  // Send across MPI channel
  else {    
    // Reset status bit
    int* status_bit = get_status_bit(direction, x, y);
    status_bit[0] = 0;
    // Add to send queue
    message msg = {x, y, direction, message_size};
    add_to_send_queue(msg, msg_buffer);
    int temp_x = x;
    int temp_y = y;
    step_in_dir(direction, &temp_x, &temp_y);
    wait_for_network(status_bit);
  }
}

char* receive_message(int x, int y, int direction){
  cpu* proc = get_cpu(x, y);
  cpu_port* port = get_port(proc, direction);
  return receive_port(port);
}

void free_network_port(int x, int y, int direction, void* buffer){
  cpu* proc = get_cpu(x, y);
  cpu_port* port = get_port(proc, direction);
  free_port(port, buffer);
}

//================================================================================
//============= Network Initialization, Start, and End ===========================
//================================================================================

void init_network(int nx, int ny){
  // Num Cpus
  nx_cpu = nx;
  ny_cpu = ny;
  int n_cpu = nx_cpu * ny_cpu;
  // Allocate and Init
  cpus = (cpu*)malloc(n_cpu * sizeof(cpu));
  for(int i=0; i<n_cpu; i++)
    init_cpu(&cpus[i]);
}

typedef struct {
  void (*thread)(int x, int y);
  int x, y;
} cpu_thread_args;

void* forwarding_thread(void* args){
  cpu_thread_args* thread_args = (cpu_thread_args*)args;
  thread_args->thread(thread_args->x, thread_args->y);
  free(thread_args);
  return 0;
}

void start_network(void (*network_thread)(int, int)){
  // For each x, y position
  for(int x=0; x<nx_cpu; x++)
    for(int y=0; y<ny_cpu; y++){
      // create thread arguments
      cpu_thread_args* args = (cpu_thread_args*)malloc(sizeof(cpu_thread_args));
      args->thread = network_thread;
      args->x = x + proc_x * nx_cpu;
      args->y = y + proc_y * ny_cpu;
      // start the thread
      start_cpu(&(cpus[x + y*nx_cpu]), &forwarding_thread, args);
    }
}

void end_network(){
  for(int i=0; i<nx_cpu*ny_cpu; i++)
    pthread_join(cpus[i].thread, NULL);
}

//================================================================================
//=============================== Main Driver ====================================
//================================================================================

// Provided Simulation
void get_simulation_size(int* num_x, int* num_y);
void network_simulation(int x, int y);
void initialize_simulation();
void finalize_simulation();

int main(int argc, char* argv[]){
  // Initialize MPI
  int rank, thread_support, n_proc;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &thread_support);
  MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // Determine Position and Size
  nx_proc = sqrt(n_proc);
  ny_proc = nx_proc;
  if(nx_proc * ny_proc != n_proc){
    printf("The specified number of processors (%d) is not a perfect square.\n", n_proc);
    exit(-1);
  }
  proc_y = rank / nx_proc;
  proc_x = rank % nx_proc;
  
  // Create Network
  get_simulation_size(&num_x, &num_y);
  init_network(num_x / nx_proc, num_y / ny_proc);

  // Initialize Simulation
  initialize_simulation();
  MPI_Barrier(MPI_COMM_WORLD);
  
  // Start MPI Thread
  init_mpi_thread();
  start_mpi_thread();

  // Start virtual network  
  start_network(&network_simulation);
  end_network();
  
  // Stop MPI Thread
  end_mpi_thread();

  // Finalize Simulation
  MPI_Barrier(MPI_COMM_WORLD);
  finalize_simulation();

  // Cleanup MPI
  MPI_Finalize();
}
