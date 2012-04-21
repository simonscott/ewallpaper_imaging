#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include "virtual_network.h"
#include "sar.h"

// Dummy row and column send operations
void send_message(int x, int y, char* message, int size){
  send_virtual_message(x * Ny + y, message, size);
}

void send_row(int ant_x, int ant_y, char* message, int size){
  for(int x = 0; x < Nx; x++)
    send_message(x, ant_y, message, size);
}

void send_col(int ant_x, int ant_y, char* message, int size){
  for(int y = 0; y < Ny; y++)
    send_message(ant_x, y, message, size);
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

void* sar_main(int MYTHREAD){
  // Get Position
  int ant_x = MYTHREAD / Ny;
  int ant_y = MYTHREAD % Ny;

  // Allocate data buffer
  int data_buffer_size = 2*sizeof(int) + Nf*sizeof(complex);
  char* data_buffer = (char*)safe_malloc(data_buffer_size, "Could not allocate data buffer.");
  complex* s = (complex*)(data_buffer + 2*sizeof(int));
  
  // Generate some fake data
  int* id = (int*)data_buffer;
  id[0] = MYTHREAD;
  for(int i=0; i<Nf; i++){
    int* words = (int*)&s[i];
    words[0] = (ant_x << 16) + ant_y;
    words[1] = i;
  }

  // Step 1: Exchange data with all processors in the same row
  // i. send data
  int* tag = (int*)data_buffer;
  tag[0] = ant_x;
  tag[1] = ant_y;
  send_row(ant_x, ant_y, (char*)tag, data_buffer_size);
  
  // ii. receive data
  for(int i=0; i<Nx; i++){
    char* msg = receive_virtual_message(MYTHREAD);
    int* msg_tag = (int*)msg;
    int source_x = msg_tag[0];
    int source_y = msg_tag[1];

    complex* msg_s = (complex*)(msg + 2*sizeof(int));
    s[source_x] = msg_s[ant_x];
    free_virtual_message(MYTHREAD, msg);
  }


    if(ant_x == 2 && ant_y == 4)
      dump_data_buffer(data_buffer);

  // Start doing data migration ...
  // printf("ant_x = %d, ant_y = %d\n", ant_x, ant_y);
}

int main(){
  start_virtual_network(Nx, Ny, &sar_main);
  return 0;
}
