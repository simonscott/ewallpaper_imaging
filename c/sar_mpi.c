#include "sar.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <pthread.h>

//================================================================================
//==================== Shared Memory Pointers ====================================
//================================================================================

complex* shared_s;                        // The data to process
complex* shared_pkt_buf;                  // The packet buffer
int N_proc_x, N_proc_y, proc_x, proc_y;   // MPI processor grid
int N_thread_x, N_thread_y;               // Num threads per core

//================================================================================
//==================== Network Functions =========================================
//================================================================================

// Sends the provided message to all processors in the same row as this processor.
// This is implemented using MPI Broadcast.
void send_row(void* message, int size)
{
}

// Sends the provided message to all processors in the same column as this processor.
// This is implemented using MPI Broadcast.
void send_col(void* message, int size)
{
}

// Receives a message from any processor in the same row as this processor.
// The message is placed in the provided buffer, and the size of the message (in bytes)
// is returned.
// This is implemented using an MPI Blocking Read.
int receive_row(void* buffer)
{
  return 0;
}

// Receives a message from any processor in the same column as this processor.
// The message is placed in the provided buffer, and the size of the message (in bytes)
// is returned.
// This is implemented using an MPI Blocking Read.
int receive_col(void* buffer)
{
  return 0;
}

//================================================================================
//==================== Helper Functions ==========================================
//================================================================================

void check_mpi_result(int status)
{
  if(status == MPI_SUCCESS)
    return;
  else {
    printf("MPI error code %d\n", status);
    exit(-1);
  }
}

//================================================================================
//==================== Main Processing Thread ====================================
//================================================================================
void* chip_thread(void* threadid)
{
  long tid = (long)(threadid);

  // Calculate thread x and y within this processor
  int thread_x = tid / N_thread_y;
  int thread_y = tid - thread_x * N_thread_y;

  // Calculate antenna x and y index in entire array
  int ant_x = proc_x * N_thread_x + thread_x;
  int ant_y = proc_y * N_thread_y + thread_y;

  // Compute pointers to own portion of shared memory
  complex* s = shared_s + tid * Nf;
  complex* pkt_buf = shared_pkt_buf + tid * (Nf + 16);

  // Actually run the SAR algorithm
  printf("%d %d\n", ant_x, ant_y);

  pthread_exit(NULL);
}

//================================================================================
//==================== The Main Function =========================================
//================================================================================

int main(int argc, char* argv[])
{
  // Declare local variables
  int N_proc, rank;
  int i, j, f, res, s_idx;
  pthread_t* threads;
  pthread_attr_t attrib;
  void* status;

  // Local arrays and buffers
  complex *recv_buf, *file_buf;              // For file ops on Node 0
  
  // MPI setup
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &N_proc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Datatype COMPLEX;
  MPI_Type_contiguous(2, MPI_FLOAT, &COMPLEX);
  MPI_Type_commit(&COMPLEX);

  // Determine position in 2D processor array
  // Assume that the processors are y-major ordered
  N_proc_x = N_proc_y = int(sqrt(N_proc));
  
  if(N_proc_x * N_proc_y != N_proc) {
    printf("Error: N_proc is not a perfect square!\n");
    return -1;
  }

  proc_x = rank / N_proc_y;
  proc_y = rank - proc_x * N_proc_y;

  // Since there are fewer physical cores than eWallpaper chips, we need
  // to use Pthreads to simulate the rest
  N_thread_x = Nx / N_proc_x;
  N_thread_y = Ny / N_proc_y;

  // Create local memory for each simulated wallpaper chip
  shared_s = (complex*)safe_malloc(Nf * N_thread_x * N_thread_y * sizeof(complex),
                  "Failed to malloc memory for s");
  shared_pkt_buf = (complex*)safe_malloc((Nf + 16) * N_thread_x * N_thread_y * sizeof(complex),
                  "Failed to malloc memory for pkt_buf");

  file_buf = (complex*)safe_malloc(Nx * Ny * Nf * sizeof(complex),
                  "Failed to malloc memory for file_buf");
  recv_buf = (complex*)safe_malloc(Nx * Ny * Nf * sizeof(complex),
                  "Failed to malloc memory for recv_buf");

  // Node 0 reads the input file and broadcasts data to all other nodes
  if(rank == 0)
    read_data(file_buf, "scene_4.dat");

  res = MPI_Bcast(file_buf, Nx * Ny * Nf, COMPLEX, 0, MPI_COMM_WORLD);
  check_mpi_result(res);

  // Extract just the data from the file that the local simulated chips require
  s_idx = 0;

  for(i = 0; i < N_thread_x; i++)
  {
    int thread_x = proc_x * N_thread_x + i;

    for(j = 0; j < N_thread_y; j++)
    {
      int thread_y = proc_y * N_thread_y + j;

      for(f = 0; f < Nf; f++)
      {
        shared_s[s_idx] = file_buf[thread_x * Ny * Nf + thread_y * Nf + f];
        s_idx++;
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // Launch the threads
  threads = (pthread_t*)safe_malloc(N_thread_x * N_thread_y * sizeof(pthread_t),
                        "Failed to malloc space for threads");
  pthread_attr_init(&attrib);
  pthread_attr_setdetachstate(&attrib, PTHREAD_CREATE_JOINABLE);

  for(i = 0; i < N_thread_x; i++) {
    for(j = 0; j < N_thread_y; j++) {
      res = pthread_create(&threads[i*N_thread_y + j], &attrib, chip_thread, (void*)(i*N_thread_y + j));

      if(res) {
        printf("Failed to create threads\n");
        exit(-1);
      }
    }
  }

  // Wait until all threads have completed
  for(i = 0; i < N_thread_x; i++) {
    for(j = 0; j < N_thread_y; j++) {
      res = pthread_join(threads[i*N_thread_y + j], &status);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // Node 0 gathers the results from all other nodes
  MPI_Gather(shared_s, N_thread_x * N_thread_y * Nf, COMPLEX, recv_buf, N_thread_x * N_thread_y * Nf, COMPLEX, 0, MPI_COMM_WORLD);

  // Node 0 writes the results to file
  if(rank == 0)
  {
    // First rearrange the data
    s_idx = 0;
    for(int mpi_proc_x = 0; mpi_proc_x < N_proc_x; mpi_proc_x++) {
      for(int pthread_x = 0; pthread_x < N_thread_x; pthread_x++) {
        for(int mpi_proc_y = 0; mpi_proc_y < N_proc_y; mpi_proc_y++) {
          for(int pthread_y = 0; pthread_y < N_thread_y; pthread_y++) {

            int offset = mpi_proc_x * N_proc_y * N_thread_x * N_thread_y * Nf +
                         mpi_proc_y * N_thread_x * N_thread_y * Nf +
                         pthread_x * N_thread_y * Nf + pthread_y * Nf;

            for(f = 0; f < Nf; f++)
            {
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
  pthread_attr_destroy(&attrib);
  free(shared_s);
  free(shared_pkt_buf);
  free(file_buf);
  free(recv_buf);
  free(threads);

  MPI_Finalize();
  pthread_exit(NULL);
}
