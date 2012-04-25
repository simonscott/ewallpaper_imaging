#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <pthread.h>
#include <math.h>
#include <string.h>
#include "sar.h"
#include "virtual_network.h"

//================================================================================
//================== Network Parameters ==========================================
//================================================================================
// N_thread_x, and N_thread_y indicate the number of virtual processors per physical
// processor. N_proc_x, and N_proc_y indicate the number of physical processors.
// proc_x, and proc_y indicate the x and y coordinate of the current physical
// processor.
int proc_x, proc_y, proc;
int N_proc_x, N_proc_y;
int N_thread_x, N_thread_y;

//================================================================================
//======================= MPI Thread =============================================
//================================================================================

int rank, N_proc;

void* mpi_thread(void* args){
  printf("mpi_thread\n");
  printf("mpi_thread: inbox_top = '%c'\n", processors[0].inbox_top[0]);
  return 0;
}

//================================================================================
//======================== Network Functions =====================================
//================================================================================

void send_message(int src_ant_x, int src_ant_y, int dest_ant_x, int dest_ant_y,
                  char* message, int size, char* send_buf){
  printf("send_message\n");
  // Encode src_x and src_y
  int header_size = 2*sizeof(int);
  int* tag = (int*)(send_buf);
  tag[0] = src_ant_x;
  tag[1] = src_ant_y;

  // Copy message
  memcpy(send_buf + header_size, message, size);

  int dest_proc_x = dest_ant_x / N_thread_x;
  int dest_proc_y = dest_ant_y / N_thread_y;
  if(dest_proc_x == proc_x && dest_proc_y == proc_y){
    int thread_x = dest_ant_x - proc_x * N_thread_x;
    int thread_y = dest_ant_y - proc_y * N_thread_y;
    int threadid = thread_x * N_thread_y + thread_y;
    printf("threadid = %d\n", threadid);
    send_virtual_message(threadid, send_buf, size + header_size);
  }
  else {
    int dest_rank = dest_proc_x * N_proc_y + dest_proc_y;
    int dest_tag = (dest_ant_x << 16) + dest_ant_y;
    printf("dest_tank = %d\n", dest_rank);
    MPI_Send(send_buf, size + header_size, MPI_BYTE, dest_rank, dest_tag, MPI_COMM_WORLD);
  }
}

char* receive_message(int threadid, int* src_x, int* src_y, int* size){
  char* msg = receive_virtual_message(threadid);
  int header_size = 2*sizeof(int);
  int* tag = (int*)(msg);
  *src_x = tag[0];
  *src_y = tag[1];
  *size = get_message_size(threadid) - header_size;
  return msg + header_size;
}

// free the message
// just forward command to virtual network.
void free_message(int threadid){
  free_virtual_message(threadid);
}

//================================================================================
//=================== Send row and Send Column ===================================
//================================================================================

void send_row(int ant_x, int ant_y, char* message, int size, char* send_buf){
  if(ant_x > 0)
    send_message(ant_x, ant_y, ant_x-1, ant_y, message, size, send_buf);  
  if(ant_x < Nx-1)
    send_message(ant_x, ant_y, ant_x+1, ant_y, message, size, send_buf);  
}

char* receive_and_forward(int ant_x, int ant_y, int* src_x, int* src_y, int* size, int threadid, char* send_buf){
  char* msg = receive_message(threadid, src_x, src_y, size);
  if(*src_x < ant_x) {
    if(ant_x < Nx - 1)
      send_message(*src_x, *src_y, ant_x+1, ant_y, msg, *size, send_buf);
  } else if(*src_x > ant_x){
    if(ant_x > 0)
      send_message(*src_x, *src_y, ant_x-1, ant_y, msg, *size, send_buf);
  } else if(*src_y < ant_y){
    if(ant_y < Ny - 1)
      send_message(*src_x, *src_y, ant_x, ant_y+1, msg, *size, send_buf);
  } else if(*src_y > ant_y){
    if(ant_y > 0)
      send_message(*src_x, *src_y, ant_x, ant_y-1, msg, *size, send_buf);
  } else {
    printf("ERROR: Unreachable Statement.\n");
  }
  return msg;
}

//================================================================================
//============================= Chip Thread ======================================
//================================================================================

void* chip_thread(int threadid){
  printf("[chip_thread] inbox_top[0] = '%c'\n", processors[0].inbox_top[0]);  
  //printf("inbox-top = '%c'\n", processors[(int)threadid].inbox_top[0]);
}

void* my_thread(void* args){
  printf("my_thread:\n");
  printf("my_thread inbox_top[0] = '%c'\n", processors[0].inbox_top[0]);
}

char* mem;

int main(int argc, char* argv[]){  
  // Initialize MPI
  // Initialize, get rank and size
  int N_proc, rank, thread_provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &thread_provided);
  MPI_Comm_size(MPI_COMM_WORLD, &N_proc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Determine Position
  N_proc_x = N_proc_y = (int)sqrt(N_proc);  
  if(N_proc_x * N_proc_y != N_proc) {
    printf("Error: N_proc is not a perfect square!\n");
    return -1;
  }
  proc_x = rank / N_proc_y;
  proc_y = rank - proc_x * N_proc_y;
  
  // Calculate thread x and y within this processor
  N_thread_x = Nx / N_proc_x;
  N_thread_y = Ny / N_proc_y;
  
  // Create the virtual network
  init_virtual_network(N_thread_x, N_thread_y);
  mem = (char*)malloc(10);
  mem[0] = 'a';
  printf("mem[0] = '%c'\n", mem[0]);
  processors[0].inbox_top[0] = 'c';
  printf("inbox_top[0] = '%c'\n", processors[0].inbox_top[0]);

  pthread_t my_thread_obj;
  pthread_create(&my_thread_obj, NULL, &my_thread, NULL);
  
  
  // Start MPI Receive Thread
  pthread_t mpi_receive_thread;
  pthread_create(&mpi_receive_thread, NULL, &mpi_thread, NULL);

  //Start the Virtual Network
  //MPI_Barrier(MPI_COMM_WORLD);
  pthread_t my_thread;
  //pthread_create(&my_thread, NULL, (void* (*)(void*))&chip_thread, NULL);
  //pthread_join(my_thread, NULL);
  start_virtual_network(&chip_thread);
  
  pthread_join(my_thread_obj, NULL);

  // Stop the MPI Thread
  //MPI_Barrier(MPI_COMM_WORLD);
  //MPI_Send(NULL, 0, MPI_BYTE, rank, 0, MPI_COMM_WORLD);
  //pthread_join(mpi_receive_thread, NULL);

  // Cleanup
  MPI_Finalize();
}
