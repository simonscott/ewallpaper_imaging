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

// send_message(src_x, src_y, dest_x, dest_y,
//              message, size,
//              send_buf)
// 1. The header size is a constant 2*sizeof(int)
// 2. Encode src_x and src_y into the tag of the send_buffer
// 3. Copy the message into the send_buffer with an offset of header_size
// 4. if dest_x and dest_y is on the current chip, then use send_virtual_message
//    i. compute the virtual threadid
// 5. otherwise, use MPI send
//    i. compute the rank of the destination processor
//    ii. pack the destination address into a single integer tag
void send_message(int src_ant_x, int src_ant_y, int dest_ant_x, int dest_ant_y,
                  char* message, int size, char* send_buf){
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
    send_virtual_message(threadid, send_buf, size + header_size);
  }
  else {
    int dest_rank = dest_proc_x * N_proc_y + dest_proc_y;
    int dest_tag = (dest_ant_x << 16) + dest_ant_y;
    MPI_Send(send_buf, size + header_size, MPI_BYTE, dest_rank, dest_tag, MPI_COMM_WORLD);
  }
}

// receive_message(virtual_threadid)
// just forward command to the virtual network.
// 1. Get the encoded message from the virtual network
// 2. Retrieve the source from the tag of the encoded message
// 3. Returns the size of the *decoded* message, and a pointer to the *decoded* message
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
    printf("BAD 1\n");
  } else if(*src_y > ant_y){
    if(ant_y > 0)
      send_message(*src_x, *src_y, ant_x, ant_y-1, msg, *size, send_buf);
    printf("BAD 2\n");
  } else {
    printf("ERROR: Unreachable Statement.\n");
  }
  return msg;
}
/*
char* receive_row(int ant_x, int ant_y, int* src_x, int* src_y, int* size, int threadid, char* send_buf){
  char* msg = receive_message(threadid, src_x, src_y, size);

  if(*src_x < ant_x && ant_x < Nx-1)
    send_message(*src_x, *src_y, ant_x+1, ant_y, msg, *size, send_buf);
  else if(*src_x > ant_x && ant_x > 0)
    send_message(*src_x, *src_y, ant_x-1, ant_y, msg, *size, send_buf);

  return msg; 
}
*/
void send_col(int ant_x, int ant_y, char* message, int size, char* send_buf){
  if(ant_y > 0)
    send_message(ant_x, ant_y, ant_x, ant_y-1, message, size, send_buf);  
  if(ant_y < Ny-1)
    send_message(ant_x, ant_y, ant_x, ant_y+1, message, size, send_buf);  
}
/*
char* receive_col(int ant_x, int ant_y, int* src_x, int* src_y, int* size, int threadid, char* send_buf){
  char* msg = receive_message(threadid, src_x, src_y, size);

  if(*src_y < ant_y && ant_y < Ny-1)
    send_message(*src_x, *src_y, ant_x, ant_y+1, msg, *size, send_buf);
  else if(*src_y > ant_y && ant_y > 0)
    send_message(*src_x, *src_y, ant_x, ant_y-1, msg, *size, send_buf);

  return msg; 
}
*/

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
    int thread_x = (dest_x - proc_x*N_thread_x);
    int thread_y = (dest_y - proc_y*N_thread_y);
    int threadid = thread_x * N_thread_y + thread_y;
    send_virtual_message(threadid, message_buffer, size);
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

  // Allocate space for send buffer
  char* send_buf = (char*)safe_malloc(MAX_MSG_SZ, "Failed to malloc memory for send buffer");

  // Start doing SAR work
  // Locate my memory
  complex* s_local_1 = shared_s + threadid * Nf;
  //complex* s_local_2 = (complex*)safe_malloc(Nf * sizeof(complex), "Failed to initialize s.");

  // Generate some fake data
  for(int i=0; i<Nf; i++){
    int* words = (int*)(&s_local_1[i]);
    words[0] = (ant_x << 16) + ant_y;
    words[1] = i;
  }
  
  // Send my local data to all processors in row
  // Every processor (x,y) sends an array, s_local_1, of size Nf, where
  // s_local_1[k] corresponds to s[x, y, k]
  send_row(ant_x, ant_y, (char*)s_local_1, Nf * sizeof(complex), send_buf);

  // Receive neighbours data
  // Receive row_s[src_x, k] from neighbours
  // where row_s[src_x, k] corresponds to s[src_x, ant_y, k]
  // Extract values from row_s where src_x = 0 to 127
  //                     v1 = s[src_x, ant_y, 2*ant_x] = row_s[src_x, 2*ant_x]
  //                     v2 = s[src_x, ant_y, 2*ant_x + 1] = row_s[src_x, 2*ant_x + 1]
  // Store values
  //                     v1 -> s_local_1[src_x]
  //                     v2 -> s_local_1[src_x+128]
  
  // Copy our own data to the correct location
  //s_local_1[ant_x] = s_local_1[2*ant_x];
  //s_local_1[ant_x + Nf/2] = s_local_1[2*ant_x + 1];
  // Receive neighbours data
  int step_1_count = 0;
  int step_2_count = 0;
  while(step_1_count < Nx - 1) {
    int src_x, src_y, size;
    complex* msg = (complex*)receive_and_forward(ant_x, ant_y, &src_x, &src_y, &size, threadid, send_buf);
    
    /* if(src_y == ant_y){ */
    /*   //Step 1 */
    /*   complex v1 = msg[2*ant_x]; */
    /*   complex v2 = msg[2*ant_x + 1]; */
    /*   s_local_1[src_x] = v1; */
    /*   s_local_1[src_x+Nf/2] = v2; */
    /*   step_1_count++; */
    /* } */
    /* else if(src_x == ant_x){ */
    /*   //Step 2 */
    /*   /\*       */
    /*   complex v1 = msg[2*ant_y]; */
    /*   complex v2 = msg[2*ant_y+1]; */
    /*   s_local_2[src_y] = v1; */
    /*   s_local_2[src_y + Nf/2] = v2; */
    /*   step_2_count++; */
    /*   *\/ */
    /*   printf("ERROR: unreachable statement (step 2)\n"); */
    /* } */
    
    free_message(threadid);
    step_1_count++;
  }

  // Debugging
  if(ant_x == 42 && ant_y == 31){
    for(int i=0; i<Nf; i++){
      int* word = (int*)(&s_local_1[i]);
      printf("[x=%d, y=%d, f=%d]\n", word[0] >> 16, word[0] & 0xFFFF, word[1]);
    }
    printf("[Step 1] Data exchange across row complete\n");
  }

  // Perform a 1D FFT twice (one for each set of frequency).

  // Send my local data to all processors in my column
  // Every processor (x,y) sends an array, s_local, of size Nf, where
  // s_local[k] corresponds to :
  //    s[k, ant_y, 2*ant_x]        if k < Nf/2,
  //    s[k-128, ant_y, 2*ant_x+1]  otherwise.
  //
  // Receive neighbours data, col_s[src_y, i], from neighbours
  // where col_s[src_y, i] corresponds to :
  //    s[i, src_y, 2*ant_x]        if i < Nf/2,
  //    s[i-128, src_y, 2*ant_x+1]  otherwise.
  // Extract values from col_s where ? = 0 to 127
  //         v1 = s[
  //
  //   s[?,?,?]

  // Send my local data to all processors in column
  /*  
  send_col(ant_x, ant_y, (char*)s_local_1, Nf * sizeof(complex), send_buf);
  
  while(step_2_count < Ny - 1) {
    int src_x, src_y, size;
    complex* msg = (complex*)receive_and_forward(ant_x, ant_y, &src_x, &src_y, &size, threadid, send_buf);
    if(src_x == ant_x){
      //Step 2
      complex v1 = msg[2*ant_y];
      complex v2 = msg[2*ant_y+1];
      s_local_2[src_y] = v1;
      s_local_2[src_y + Nf/2] = v2;
      step_2_count++;
    }
    else if(src_y == ant_y){
      //Step 3
      printf("ERROR: Unreachable Statement Step 3\n");
    }
    free_message(threadid);
  }
  */
  /*
  // Copy our own data to the correct location
  s_local_1[ant_y] = s_local_1[2*ant_y];
  s_local_1[ant_y + Nf/2] = s_local_1[2*ant_y + 1];

  // Receive neighbours data
  for(int i=0; i<Ny-1; i++) {
    int src_x, src_y, size;
    complex* msg = (complex*)receive_col(ant_x, ant_y, &src_x, &src_y, &size, threadid, send_buf);
    complex v1 = msg[2*ant_y];
    complex v2 = msg[2*ant_y+1];
    s_local_1[src_y] = v1;
    s_local_1[src_y + Nf/2] = v2;
    free_message(threadid);
  }
  

  // Debugging
  if(ant_x == 42 && ant_y == 31){
    for(int i=0; i<Nf; i++){
      int* word = (int*)(&s_local_2[i]);
      printf("[x=%d, y=%d, f=%d]\n", word[0] >> 16, word[0] & 0xFFFF, word[1]);
    }
    printf("[Step 2] Data exchange across column complete\n");
  }
  */
  // Finish thread
  //free(s_local_2);
  free(send_buf);
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

  // Create virtual network
  init_virtual_network(N_thread_x, N_thread_y);

  // Start the MPI Thread
  pthread_t mpi_receive_thread;
  pthread_create(&mpi_receive_thread, NULL, &mpi_thread, NULL);

  // Start the virtual network
  start_virtual_network(chip_thread);

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
