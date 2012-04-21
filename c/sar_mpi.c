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
//===================== Network Functions ========================================
//================================================================================

// send_message(dest_x, dest_y, message, size)
// 1. if dest_x and dest_y is on the current chip, then use send_virtual_message
//    i. compute the virtual threadid
// 2. otherwise, use MPI send
//    i. compute the rank of the destination processor
//    ii. pack the destination address into a single integer tag
void send_message(int dest_ant_x, int dest_ant_y, char* message, int size){
  if ((dest_ant_x / N_thread_x == proc_x) && (dest_ant_y / N_thread_y == proc_y)){
    int threadid = (dest_ant_x - proc_x*N_thread_x) * N_thread_y + (dest_ant_y - proc_y * N_thread_y);
    send_virtual_message(threadid, message, size);
  }
  else {
    int dest_proc_x = dest_ant_x / N_thread_x;
    int dest_proc_y = dest_ant_y / N_thread_y;
    int dest_rank = dest_proc_x * N_proc_y + dest_proc_y;
    int tag = (dest_ant_x << 16) + dest_ant_y;
    MPI_Send(message, size, MPI_BYTE, dest_rank, tag, MPI_COMM_WORLD);
  }
}

// receive_message(virtual_threadid)
// just forward command to the virtual network.
char* receive_message(int threadid){
  return receive_virtual_message(threadid);
}

// free the message
// just forward command to virtual network.
void free_message(int threadid, char* message){
  free_virtual_message(threadid, message);
}

//================================================================================
//=================== Send row and Send Column ===================================
//================================================================================

void send_row(int ant_x, int ant_y, char* message, int size){
  for(int x = 0; x < Nx; x++)
    send_message(x, ant_y, message, size);
}

void send_col(int ant_x, int ant_y, char* message, int size){
  for(int y = 0; y < Ny; y++)
    send_message(ant_x, y, message, size);
}


//================================================================================
//======================= MPI Thread =============================================
//================================================================================

// mpi_thread
// continually receives messages and forwards messages to the appropriate
// virtual processor.
void* mpi_thread(void* args){
  char* message_buffer = (char*)malloc(message_memory);
  while(1){
    MPI_Status status;
    MPI_Recv(message_buffer, message_memory, MPI_BYTE,
             MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    pthread_testcancel();
    //Compute destination
    int dest_x = status.MPI_TAG >> 16;
    int dest_y = status.MPI_TAG & 0xFFFF;
    int size; MPI_Get_count(&status, MPI_BYTE, &size);
    send_message(dest_x, dest_y, message_buffer, size);
  }
}

//================================================================================
//==================== Main Processing Thread ====================================
//================================================================================

complex* shared_s;

// chip_thread(virtual_threadid)
// 1. Determine thread_x, and thread_y, the coordinates of the thread within this
//    virtual network.
// 2. Determine ant_x, and ant_y, the absolute coordinates on the whole network.

void* chip_thread(int threadid)
{
  // Calculate thread x and y within this processor
  int thread_x = threadid / N_thread_y;
  int thread_y = threadid - thread_x * N_thread_y;

  // Calculate antenna x and y index in entire array
  int ant_x = proc_x * N_thread_x + thread_x;
  int ant_y = proc_y * N_thread_y + thread_y;

  // Start doing SAR work
  // Locate my memory
  complex* s = shared_s + threadid * Nf;

  s[0].real = ant_x;
  s[0].imag = ant_y;
  send_message((ant_x + 1)%Nx, ant_y, (char*)s, sizeof(complex));

  complex* msg = (complex*)receive_message(threadid);
  s[0] = msg[0];
  free_message(threadid, (char*)msg);
  printf("antenna (%d,%d) received (%f,%f)\n", ant_x, ant_y, s[0].real, s[0].imag);

  // Generate some fake data
  int* id = (int*)data_buffer;
  id[0] = MYTHREAD;
  for(int i=0; i<Nf; i++){
    int* words = (int*)&s[i];
    words[0] = (ant_x << 16) + ant_y;
    words[1] = i;
  }

  // Finish thread
  pthread_exit(NULL);
}

//================================================================================
//====================== Main Driver =============================================
//================================================================================

// Ensure successful MPI broadcast.
void check_mpi_result(int status){
  if(status == MPI_SUCCESS)
    return;
  else {
    printf("MPI error code %d\n", status);
    exit(-1);
  }
}

// Initialize MPI
//   1. Initialize MPI, get rank and size
//   2. Declare COMPLEX type
// Determine position
//   1. Determine N_proc_x and N_proc_y
//   2. Error if we were not given a perfect square number of processors
//   3. Calculate position:
//      proc_x = rank / N_proc_y
//      proc_y = rank - proc_x * N_proc_y
// Create virtual network
//   1. Determine how many virtual processors to create.
//   2. Start the network with entry chip_thread.
int main(int argc, char* argv[]){  
  // Initialize MPI
  // Initialize, get rank and size
  int N_proc, rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &N_proc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // Declare COMPLEX type
  MPI_Datatype COMPLEX;
  MPI_Type_contiguous(2, MPI_FLOAT, &COMPLEX);
  MPI_Type_commit(&COMPLEX);
  
  // Determine Position
  N_proc_x = N_proc_y = (int)sqrt(N_proc);  
  if(N_proc_x * N_proc_y != N_proc) {
    printf("Error: N_proc is not a perfect square!\n");
    return -1;
  }
  proc_x = rank / N_proc_y;
  proc_y = rank - proc_x * N_proc_y;
  
  // Calculate thread x and y within this processor, and create the virtual network
  N_thread_x = Nx / N_proc_x;
  N_thread_y = Ny / N_proc_y;

  // Create local memory for each simulated wallpaper chip
  shared_s = (complex*)safe_malloc(Nf * N_thread_x * N_thread_y * sizeof(complex),
                  "Failed to malloc memory for s");
  complex* file_buf = (complex*)safe_malloc(Nx * Ny * Nf * sizeof(complex),
                  "Failed to malloc memory for file_buf");
  complex* recv_buf = (complex*)safe_malloc(Nx * Ny * Nf * sizeof(complex),
                  "Failed to malloc memory for recv_buf");

  // Node 0 reads the input file and broadcasts data to all other nodes
  if(rank == 0)
    read_data(file_buf, "scene_4.dat");
  int res = MPI_Bcast(file_buf, Nx * Ny * Nf, COMPLEX, 0, MPI_COMM_WORLD);
  check_mpi_result(res);

  // Extract just the data from the file that the local simulated chips require
  int s_idx = 0;
  for(int i = 0; i < N_thread_x; i++) {
    int thread_x = proc_x * N_thread_x + i;
    for(int j = 0; j < N_thread_y; j++) {
      int thread_y = proc_y * N_thread_y + j;
      for(int f = 0; f < Nf; f++) {
        shared_s[s_idx] = file_buf[thread_x * Ny * Nf + thread_y * Nf + f];
        s_idx++;
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // Start the MPI Thread
  pthread_t mpi_receive_thread;
  pthread_create(&mpi_receive_thread, NULL, &mpi_thread, NULL);

  // Create and start virtual network
  start_virtual_network(N_thread_x, N_thread_y, chip_thread);

  // Stop the MPI Thread  
  pthread_cancel(mpi_receive_thread);
  MPI_Send(NULL, 0, MPI_BYTE, rank, 0, MPI_COMM_WORLD);
  pthread_join(mpi_receive_thread, NULL);

  // Wait for simulation to finish and gather results
  MPI_Barrier(MPI_COMM_WORLD);
  // Node 0 gathers the results from all other nodes
  MPI_Gather(shared_s, N_thread_x * N_thread_y * Nf, COMPLEX, recv_buf,
             N_thread_x * N_thread_y * Nf, COMPLEX, 0, MPI_COMM_WORLD);  
  // Node 0 writes the results to file
  if(rank == 0) {
    // First rearrange the data
    s_idx = 0;
    for(int mpi_proc_x = 0; mpi_proc_x < N_proc_x; mpi_proc_x++) {
      for(int pthread_x = 0; pthread_x < N_thread_x; pthread_x++) {
        for(int mpi_proc_y = 0; mpi_proc_y < N_proc_y; mpi_proc_y++) {
          for(int pthread_y = 0; pthread_y < N_thread_y; pthread_y++) {
            int offset = mpi_proc_x * N_proc_y * N_thread_x * N_thread_y * Nf +
                         mpi_proc_y * N_thread_x * N_thread_y * Nf +
                         pthread_x * N_thread_y * Nf + pthread_y * Nf;
            for(int f = 0; f < Nf; f++) {
              file_buf[s_idx] = recv_buf[offset + f];
              s_idx++;
            }
          }
        }
      }
    }

    // Then write to file
    write_data(file_buf, "scene_4_mpi.out");
  }

  // Free allocated memory
  free(shared_s);

  // Cleanup MPI
  MPI_Finalize();
  pthread_exit(NULL);
}
