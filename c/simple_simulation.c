#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "network.h"
#include "virtual_cpu.h"
#include "sar.h"

//================================================================================
//============================ Higher Level Network Functions ====================
//================================================================================

void send_row(int x, int y, void* data, int size){
  if(x > 0)
    send_message(x, y, left, data, size);
  if(x < Nx - 1)
    send_message(x, y, right, data, size);
}

void send_col(int x, int y, void* data, int size){
  if(y > 0)
    send_message(x, y, down, data, size);
  if(y < Ny - 1)
    send_message(x, y, up, data, size);
}

typedef struct {
  int x;
  int y;
  int num_less;
  int num_more;
  int receive_dir;
  char* buffer;
  int size;
} receive_status;

receive_status start_receive_row(int x, int y, char* msg_buffer, int msg_size){
  receive_status status;
  status.x = x;
  status.y = y;
  status.num_less = x;
  status.num_more = Nx - (x + 1);
  status.receive_dir = right;
  status.buffer = msg_buffer;
  status.size = msg_size;
  return status;
}

receive_status start_receive_col(int x, int y, char* msg_buffer, int msg_size){
  receive_status status;
  status.x = x;
  status.y = y;
  status.num_less = y;
  status.num_more = Ny - (y + 1);
  status.receive_dir = up;
  status.buffer = msg_buffer;
  status.size = msg_size;
  return status;
}

char* receive_line(receive_status* status, int less_dir,
                   int coordinate, int coordinate_max){
  int more_dir = opposite_dir(less_dir);
  
  if(status->num_less == 0 ||
     (status->num_more > 0 && status->receive_dir == more_dir)){
      status->receive_dir = less_dir;
      // Receive and copy message into buffer
      //printf("(%d,%d) receiving %d\n", status->x, status->y, more_dir); //DEBUG
      char* msg = receive_message(status->x, status->y, more_dir);
      // Forward message
      memcpy(status->buffer, (char*)msg, status->size);
      free_network_port(status->x, status->y, more_dir, msg);
      if(coordinate > 0){
        //printf("(%d,%d) sending %d\n", status->x, status->y, less_dir); //DEBUG
        send_message(status->x, status->y, less_dir, status->buffer, status->size);
      }
      // Decrement counter
      status->num_more--;
      // Process message
      return status->buffer;
  }
  else if(status->num_more == 0 ||
          (status->num_less > 0 && status->receive_dir == less_dir)){
    status->receive_dir = more_dir;
    // Receive and copy message into buffer
    //printf("(%d,%d) receiving %d\n", status->x, status->y, less_dir); //DEBUG
    char* msg = receive_message(status->x, status->y, less_dir);
    // Forward message
    memcpy(status->buffer, (char*)msg, status->size);
    free_network_port(status->x, status->y, less_dir, msg);
    if(coordinate < coordinate_max - 1){
      //printf("(%d,%d) sending %d\n", status->x, status->y, more_dir); //DEBUG
      send_message(status->x, status->y, more_dir, status->buffer, status->size);
    }
    // Decrement counter
    status->num_less--;
    // Process message
    return status->buffer;
  }
  else {
    printf("No more messages!\n");
    exit(-1);
  }
}

char* receive_row(receive_status* status){
  return receive_line(status, left, status->x, Nx);
}

char* receive_col(receive_status* status){
  return receive_line(status, down, status->y, Ny);
}

//================================================================================
//====================== Simulation ==============================================
//================================================================================
void get_simulation_size(int* x, int* y){
  x[0] = Nx;
  y[0] = Ny;
}

void dump_data_buffer(char* data_buffer){
  int* id = (int*)data_buffer;
  complex* s = (complex*)(data_buffer + 2*sizeof(int));

  printf("id = %d, %d\n", id[0], id[1]);
  for(int i=0; i<Nf; i++){
    unsigned int* words = (unsigned int*)&s[i];
    printf("[%d, %d, %d]\n", words[0] >> 16, words[0] & 0xFFFF, words[1]);
  }
}

// Needed MPI_Network functions
extern int proc_x, proc_y, nx_cpus, ny_cpus, nx_proc, ny_proc;
extern void to_local_coords(int* x, int* y);


complex* file_buffer;
void initialize_simulation(){
  file_buffer = (complex*)safe_malloc(Nx * Ny * Nf * sizeof(complex),
                                      "Failed to allocate memory for file.\n");
  
  // Read in data and broadcast to everyone
  if(proc_x == 0 && proc_y == 0)
    read_data(file_buffer, "scene_4.dat");
  MPI_Bcast(file_buffer, Nx * Ny * Nf * sizeof(complex), MPI_BYTE, 0, MPI_COMM_WORLD);
}

void finalize_simulation(){
  // Gather everyone's data in processor 0
  complex* receive_buffer = (complex*)safe_malloc(Nx * Ny * Nf * sizeof(complex),
                                                  "Failed to allocate receive buffer.\n");
  MPI_Gather(file_buffer, nx_cpus * ny_cpus * Nf * sizeof(complex), MPI_BYTE,
             receive_buffer, nx_cpus * ny_cpus * Nf * sizeof(complex), MPI_BYTE,
             0, MPI_COMM_WORLD);

  // Node 0 writes data to file
  if(proc_x == 0 && proc_y == 0){
    // First rearrange the data
    s_idx = 0;
    for(int mpi_proc_x = 0; mpi_proc_x < nx_proc; mpi_proc_x++) {
      for(int pthread_x = 0; pthread_x < nx_cpus; pthread_x++) {
        for(int mpi_proc_y = 0; mpi_proc_y < ny_proc; mpi_proc_y++) {
          for(int pthread_y = 0; pthread_y < ny_cpus; pthread_y++) {
            int offset = mpi_proc_x * ny_proc * nx_cpus * ny_cpus * Nf +
                         mpi_proc_y * nx_cpus * ny_cpus * Nf +
                         pthread_x * ny_cpus * Nf + pthread_y * Nf;
            for(int f = 0; f < Nf; f++) {
              file_buffer[s_idx] = receive_buffer[offset + f];
              s_idx++;
            }
          }
        }
      }
    }

    // Then write to file
    write_data(file_buffer, "scene_4_mpi.out");
  }
}

void network_simulation(int x, int y){
  //Initialize Buffers
  int msg_size = Nf * sizeof(complex) + 2 * sizeof(int);
  char* forwarding_buffer = (char*)malloc(msg_size);
  char* up_buffer = (char*)malloc(msg_size);
  char* down_buffer = (char*)malloc(msg_size);
  char* left_buffer = (char*)malloc(msg_size);
  char* right_buffer = (char*)malloc(msg_size);
  free_network_port(x, y, up, up_buffer);
  free_network_port(x, y, down, down_buffer);
  free_network_port(x, y, left, left_buffer);
  free_network_port(x, y, right, right_buffer);

  //Initialize Local Data
  char* data = (char*)malloc(Nf * sizeof(complex) + 2*sizeof(int));
  int* tag = (int*)data;
  tag[0] = x;
  tag[1] = y;

  //Generate some fake data
  /*
  complex* s = (complex*)(data + 2*sizeof(int));
  for(int i=0; i<Nf; i++){
    int* words = (int*)(&s[i]);
    words[0] = (x << 16) + y;
    words[1] = i;
  }
  */

  //Extract data from file buffer
  complex* s = (complex*)(data + 2*sizeof(int));
  for(int i=0; i<Nf; i++)
    s[i] = file_buffer[x*Ny*Nf + y*Nf + i];

  //--------------------------------------------------------------------------------
  //---------------------- Step 1 --------------------------------------------------

  //Send row
  send_row(x, y, data, msg_size);

  //Copy our own data
  s[x] = s[2*x];
  s[x + Nf/2] = s[2*x + 1];

  //Receive row
  receive_status status = start_receive_row(x, y, forwarding_buffer, msg_size);
  for(int i=0; i<Nx-1; i++){
    char* msg = receive_row(&status);
    int* tag = (int*)msg;
    complex* s2 = (complex*)(msg + 2*sizeof(int));

    int src_x = tag[0];
    int src_y = tag[1];
    complex v1 = s2[2*x];
    complex v2 = s2[2*x + 1];
    s[src_x] = v1;
    s[src_x + Nf/2] = v2;
  }

  if(x == 1 && y == 3)
    dump_data_buffer(data);

  //--------------------------------------------------------------------------------
  //------------------------ Step 2 ------------------------------------------------

  // Send out column
  send_col(x, y, data, msg_size);

  // Receive column
  status = start_receive_col(x, y, forwarding_buffer, msg_size);
  for(int i=0; i<Ny-1; i++){
    char* msg = receive_col(&status);
    int* tag = (int*)msg;
    complex* s2 = (complex*)(msg + 2*sizeof(int));

    int src_x = tag[0];
    int src_y = tag[1];
    complex v1 = s2[y];
    complex v2 = s2[y + Nf/2];
    s[src_y] = v1;
    s[src_y + Nf/2] = v2;
  }

  if(x == 1 && y == 3)
    dump_data_buffer(data);

  //--------------------------------------------------------------------------------
  //--------------------------- Step 3 ---------------------------------------------

  // Send row
  send_row(x, y, data, msg_size);

  // Copy our own data
  complex v1 = s[x];
  complex v2 = s[x + Nf/2];
  s[2*x] = v1;
  s[2*x+1] = v2;

  // Receive neighbours data
  status = start_receive_row(x, y, forwarding_buffer, msg_size);
  for(int i=0; i<Nx-1; i++){
    char* msg = receive_row(&status);
    int* tag = (int*)msg;
    complex* s2 = (complex*)(msg + 2*sizeof(int));

    int src_x = tag[0];
    int src_y = tag[1];
    v1 = s2[x];
    v2 = s2[x + Nf/2];
    s[2*src_x] = v1;
    s[2*src_x+1] = v2;
  }

  if(x == 1 && y == 3)
    dump_data_buffer(data);

  //--------------------------------------------------------------------------------
  //------------------------ Step 4 ------------------------------------------------
  
  // Send out column
  send_col(x, y, data, msg_size);

  // Copy our own data
  v1 = s[2*y];
  v2 = s[2*y+1];
  s[y] = v1;
  s[y+Nf/2] = v2;

  // Receive column
  status = start_receive_col(x, y, forwarding_buffer, msg_size);
  for(int i=0; i<Ny-1; i++){
    char* msg = receive_col(&status);
    int* tag = (int*)msg;
    complex* s2 = (complex*)(msg + 2*sizeof(int));

    int src_x = tag[0];
    int src_y = tag[1];
    complex v1 = s2[2*y];
    complex v2 = s2[2*y+1];
    s[src_y] = v1;
    s[src_y + Nf/2] = v2;
  }

  if(x == 1 && y == 3)
    dump_data_buffer(data);

  //--------------------------------------------------------------------------------
  //--------------------------- Step 5 ---------------------------------------------

  // Send row
  send_row(x, y, data, msg_size);

  // Receive neighbours data
  status = start_receive_row(x, y, forwarding_buffer, msg_size);
  for(int i=0; i<Nx-1; i++){
    char* msg = receive_row(&status);
    int* tag = (int*)msg;
    complex* s2 = (complex*)(msg + 2*sizeof(int));

    int src_x = tag[0];
    int src_y = tag[1];
    v1 = s2[x];
    v2 = s2[x + Nf/2];
    s[src_x] = v1;
    s[src_x+Nf/2] = v2;
  }

  if(x == 1 && y == 3)
    dump_data_buffer(data);

  //--------------------------------------------------------------------------------
  //------------------------ Step 6 ------------------------------------------------
  
  // Send out column
  send_col(x, y, data, msg_size);

  // Copy our own data
  v1 = s[y];
  v2 = s[y + Nf/2];
  s[2*y] = v1;
  s[2*y + 1] = v2;

  // Receive column
  status = start_receive_col(x, y, forwarding_buffer, msg_size);
  for(int i=0; i<Ny-1; i++){
    char* msg = receive_col(&status);
    int* tag = (int*)msg;
    complex* s2 = (complex*)(msg + 2*sizeof(int));

    int src_x = tag[0];
    int src_y = tag[1];
    complex v1 = s2[y];
    complex v2 = s2[y+Nf/2];
    s[2*src_y] = v1;
    s[2*src_y + 1] = v2;
  }

  if(x == 1 && y == 3)
    dump_data_buffer(data);

  //--------------------------------------------------------------------------------
  //----------------------- Write data to file buffer ------------------------------

  // Get local coordinates
  int local_x = x;
  int local_y = y;
  to_local_coords(&local_x, &local_y);

  for(int i=0; i<Nf; i++)
    filebuf[local_x*Ny*Nf + local_y*Nf + i] = s[i];
}
