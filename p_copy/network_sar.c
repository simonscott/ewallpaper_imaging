#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <string.h>
#include "virtual_network.h"
#include "sar.h"

//================================================================================
//======================== Nework Functions ======================================
//================================================================================

// Dummy row and column send operations
void send_message(int src_x, int src_y, int dest_x, int dest_y,
                  char* message, int size, char* send_buf){
  //Encode
  int header_size = 2*sizeof(int);
  int* tag = (int*)send_buf;
  tag[0] = src_x;
  tag[1] = src_y;
  //Copy Message
  memcpy(send_buf + header_size, message, size);
  //Send
  send_virtual_message(dest_x * Ny + dest_y, send_buf, size + header_size);
}

char* receive_message(int threadid, int* src_x, int* src_y, int* size){
  char* msg = receive_virtual_message(threadid);
  //Decode tag
  int header_size = 2*sizeof(int);
  int* tag = (int*)msg;
  *src_x = tag[0];
  *src_y = tag[1];
  *size = get_message_size(threadid) - header_size;
  return msg + header_size;
}

void free_message(int threadid){
  free_virtual_message(threadid);
}

//================================================================================
//========================= Send/Receive Row/Col =================================
//================================================================================

void send_row(int src_x, int src_y, char* message, int size, char* send_buf){
  if(src_x > 0)
    send_message(src_x, src_y, src_x - 1, src_y, message, size, send_buf);
  if(src_x < Nx-1)
    send_message(src_x, src_y, src_x + 1, src_y, message, size, send_buf);
}

void send_col(int src_x, int src_y, char* message, int size, char* send_buf){
  printf("[send_col] Not yet implemented.\n");
  exit(-1);  
}

char* receive_and_forward(int ant_x, int ant_y, int* src_x, int* src_y, int* size,
                          int threadid, char* send_buf){
  char* msg = receive_message(threadid, src_x, src_y, size);
  //Forward Right
  if(*src_x < ant_x){
    if(ant_x < Nx - 1)
      send_message(*src_x, *src_y, ant_x + 1, ant_y, msg, *size, send_buf);
  }
  //Forward Left
  else if(*src_x > ant_x){
    if(ant_x > 0)
      send_message(*src_x, *src_y, ant_x - 1, ant_y, msg, *size, send_buf);
  }
  //Forward Up
  else if(*src_y < ant_y){
    if(ant_y < Ny - 1)
      send_message(*src_x, *src_y, ant_x, ant_y+1, msg, *size, send_buf);
  }
  //Forward Down
  else if(*src_y > ant_y){
    if(ant_y > 0)
      send_message(*src_x, *src_y, ant_x, ant_y-1, msg, *size, send_buf);
  }
  //Diagonal
  else {
    printf("[receive_and_forward] Diagonal message send not supported.\n");
    exit(-1);
  }
  return msg;
}

//================================================================================
//================================================================================
//================================================================================

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

  // Allocate space for send buffer
  char* send_buf = (char*)safe_malloc(MAX_MSG_SZ, "Failure to allocate send buffer.");

  // Generate some fake data
  complex* s_local_1 = (complex*)safe_malloc(Nf * sizeof(complex), "Failed to initialize s.");
  for(int i=0; i<Nf; i++){
    int* words = (int*)(&s_local_1[i]);
    words[0] = (ant_x << 16) + ant_y;
    words[1] = i;
  }

  // Send my local data
  send_row(ant_x, ant_y, (char*)s_local_1, Nf * sizeof(complex), send_buf);

  // Receive
  for(int i=0; i<Nx-1; i++){
    int src_x, src_y, size;
    char* msg = receive_and_forward(ant_x, ant_y, &src_x, &src_y, &size, MYTHREAD, send_buf);
    printf("(%d,%d) received %d bytes from (%d,%d)\n", ant_x, ant_y, size, src_x, src_y);
    free_message(MYTHREAD);
  }
}

int main(){
  init_virtual_network(Nx, Ny);
  start_virtual_network(&sar_main);
  return 0;
}
