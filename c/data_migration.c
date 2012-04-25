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
  if(src_y > 0)
    send_message(src_x, src_y, src_x, src_y-1, message, size, send_buf);
  if(src_y < Ny-1)
    send_message(src_x, src_y, src_x, src_y+1, message, size, send_buf);
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

void* sar_main(int threadid){
  // Get Position
  int ant_x = threadid / Ny;
  int ant_y = threadid % Ny;

  // Allocate space for send buffer
  char* send_buf = (char*)safe_malloc(MAX_MSG_SZ, "Failure to allocate send buffer.");

  // Allocate space 
  complex* s_local_1 = (complex*)safe_malloc(Nf * sizeof(complex), "Failed to initialize s.");
  complex* s_local_2 = (complex*)safe_malloc(Nf * sizeof(complex), "Failed to initialize s.");

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

  // Copy our own data to the correct location
  s_local_1[ant_x] = s_local_1[2*ant_x];
  s_local_1[ant_x + Nf/2] = s_local_1[2*ant_x + 1];
  
  // Receive neighbours data
  int step_1_count = 0;
  int step_2_count = 0;
  int step_3_count = 0;
  int step_4_count = 0;
  int step_5_count = 0;
  int step_6_count = 0;
  
  while(step_1_count < Nx - 1) {
    int src_x, src_y, size;
    complex* msg = (complex*)receive_and_forward(ant_x, ant_y, &src_x, &src_y, &size, threadid, send_buf);

    if(src_y == ant_y){
      //Step 1
      complex v1 = msg[2*ant_x];
      complex v2 = msg[2*ant_x + 1];
      s_local_1[src_x] = v1;
      s_local_1[src_x+Nf/2] = v2;
      step_1_count++;
    }
    else if(src_x == ant_x){
      //Step 2
      complex v1 = msg[ant_y];
      complex v2 = msg[ant_y+Nf/2];
      s_local_2[src_y] = v1;
      s_local_2[src_y + Nf/2] = v2;
      step_2_count++;
    }

    free_message(threadid);
  }

  // Do 2 1D FFTs (one for each frequency band) each of Nx points

  // Print out some values
  if(ant_x == 7 && ant_y == 20){
    printf("Finished Data for (7,20):\n");
    for(int i=0; i<Nf; i++){
      int* word = (int*)(&s_local_1[i]);
      printf("[x=%d, y=%d, f=%d]\n", word[0] >> 16, word[0] & 0xFFFF, word[1]);
    }
    printf("[Step 1] Data exchange across row complete\n");
  }

  // Send my local data to all processors in column
  send_col(ant_x, ant_y, (char*)s_local_1, Nf * sizeof(complex), send_buf);

  // Copy out own data
  s_local_2[ant_y] = s_local_1[ant_y];
  s_local_2[ant_y + Nf/2] = s_local_1[ant_y + Nf/2];

  // Receive neighbours data
  while(step_2_count < Ny - 1) {
    int src_x, src_y, size;
    complex* msg = (complex*)receive_and_forward(ant_x, ant_y, &src_x, &src_y, &size, threadid, send_buf);
    if(src_x == ant_x){
      //Step 2
      complex v1 = msg[ant_y];
      complex v2 = msg[ant_y+Nf/2];
      s_local_2[src_y] = v1;
      s_local_2[src_y + Nf/2] = v2;
      step_2_count++;
    }
    else if(src_y == ant_y){
      //Step 3
      complex v1 = msg[ant_x];
      complex v2 = msg[ant_x + Nf/2];
      s_local_1[2*src_x] = v1;
      s_local_1[2*src_x+1] = v2;
      step_3_count++;
    }
    free_message(threadid);
  }

  // Print out some values
  if(ant_x == 7 && ant_y == 20){
    for(int i=0; i<Nf; i++){
      int* word = (int*)(&s_local_2[i]);
      printf("[x=%d, y=%d, f=%d]\n", word[0] >> 16, word[0] & 0xFFFF, word[1]);
    }
    printf("[Step 2] Data exchange across columns complete\n");
  }

  // Do 2 1D FFT (one for each frequency band) each of Ny points
  // Multiply each element by complex exponential

  
  // Send my local data to all processors in row
  send_row(ant_x, ant_y, (char*)s_local_2, Nf * sizeof(complex), send_buf);

  // Copy out own data
  s_local_1[2*ant_x] = s_local_2[ant_x];
  s_local_1[2*ant_x+1] = s_local_2[ant_x + Nf/2];

  // Receive neighbours data
  while(step_3_count < Nx - 1) {
    int src_x, src_y, size;
    complex* msg = (complex*)receive_and_forward(ant_x, ant_y, &src_x, &src_y, &size, threadid, send_buf);
    if(src_y == ant_y){
      //Step 3
      complex v1 = msg[ant_x];
      complex v2 = msg[ant_x + Nf/2];
      s_local_1[2*src_x] = v1;
      s_local_1[2*src_x+1] = v2;
      step_3_count++;
    }
    else if(src_x == ant_x){
      //Step 4
      complex v1 = msg[2*ant_y];
      complex v2 = msg[2*ant_y+1];
      s_local_2[src_y] = v1;
      s_local_2[src_y + Nf/2] = v2;
      step_4_count++;
    }
    free_message(threadid);
  }  

  // Print out some values
  if(ant_x == 7 && ant_y == 20){
    for(int i=0; i<Nf; i++){
      int* word = (int*)(&s_local_1[i]);
      printf("[x=%d, y=%d, f=%d]\n", word[0] >> 16, word[0] & 0xFFFF, word[1]);
    }
    printf("[Step 3] Data exchange across rows complete\n");
  }


  // Do linear interpolation over the entire array of Nf points
  // Do a 1D IFFT over entire array of Nf points


  
  send_col(ant_x, ant_y, (char*)s_local_1, Nf * sizeof(complex), send_buf);

  // Copy out own data
  s_local_2[ant_y] = s_local_1[2*ant_y];
  s_local_2[ant_y + Nf/2] = s_local_1[2*ant_y+1];

  // Receive neighbours data
  while(step_4_count < Ny - 1) {
    int src_x, src_y, size;
    complex* msg = (complex*)receive_and_forward(ant_x, ant_y, &src_x, &src_y, &size, threadid, send_buf);
    if(src_x == ant_x){
      //Step 4
      complex v1 = msg[2*ant_y];
      complex v2 = msg[2*ant_y+1];
      s_local_2[src_y] = v1;
      s_local_2[src_y + Nf/2] = v2;
      step_4_count++;
    }
    else if(src_y == ant_y){
      //Step 5
      complex v1 = msg[ant_x];
      complex v2 = msg[ant_x + Nf/2];
      s_local_1[src_x] = v1;
      s_local_1[src_x+Nf/2] = v2;
      step_5_count++;
    }
    free_message(threadid);
  }

  // Print out some values
  if(ant_x == 7 && ant_y == 20){
    for(int i=0; i<Nf; i++){
      int* word = (int*)(&s_local_2[i]);
      printf("[x=%d, y=%d, f=%d]\n", word[0] >> 16, word[0] & 0xFFFF, word[1]);
    }
    printf("[Step 4] Data exchange across columns complete\n");
  }


  // Do 2 1D FFTs, each of Nx points

  send_row(ant_x, ant_y, (char*)s_local_2, Nf * sizeof(complex), send_buf);

  // Copy out own data
  s_local_1[ant_x] = s_local_2[ant_x];
  s_local_1[ant_x + Nf/2] = s_local_2[ant_x + Nf/2];

  // Receive neighbours data
  while(step_5_count < Nx - 1) {
    int src_x, src_y, size;
    complex* msg = (complex*)receive_and_forward(ant_x, ant_y, &src_x, &src_y, &size, threadid, send_buf);
    if(src_y == ant_y){
      //Step 5
      complex v1 = msg[ant_x];
      complex v2 = msg[ant_x + Nf/2];
      s_local_1[src_x] = v1;
      s_local_1[src_x+Nf/2] = v2;
      step_5_count++;
    }
    else if(src_x == ant_x){
      //Step 6
      complex v1 = msg[ant_y];
      complex v2 = msg[ant_y+Nf/2];
      s_local_2[2*src_y] = v1;
      s_local_2[2*src_y + 1] = v2;
      step_6_count++;
    }
    free_message(threadid);
  }  

  // Print out some values
  if(ant_x == 7 && ant_y == 20){
    for(int i=0; i<Nf; i++){
      int* word = (int*)(&s_local_1[i]);
      printf("[x=%d, y=%d, f=%d]\n", word[0] >> 16, word[0] & 0xFFFF, word[1]);
    }
    printf("[Step 5] Data exchange across rows complete\n");
  }

  // Do 2 1D FFts, each of Ny points

  send_col(ant_x, ant_y, (char*)s_local_1, Nf * sizeof(complex), send_buf);

  // Copy out own data
  s_local_2[2*ant_y] = s_local_1[ant_y];
  s_local_2[2*ant_y + 1] = s_local_1[ant_y + Nf/2];

  // Receive neighbours data
  while(step_6_count < Ny - 1) {
    int src_x, src_y, size;
    complex* msg = (complex*)receive_and_forward(ant_x, ant_y, &src_x, &src_y, &size, threadid, send_buf);
    if(src_x == ant_x){
      //Step 6
      complex v1 = msg[ant_y];
      complex v2 = msg[ant_y+Nf/2];
      s_local_2[2*src_y] = v1;
      s_local_2[2*src_y + 1] = v2;
      step_6_count++;
    }
    else if(src_y == ant_y){
      //Step 7
      printf("Unreachable Statement\n");
    }
    free_message(threadid);
  }

  // Print out some values
  if(ant_x == 7 && ant_y == 20){
    for(int i=0; i<Nf; i++){
      int* word = (int*)(&s_local_2[i]);
      printf("[x=%d, y=%d, f=%d]\n", word[0] >> 16, word[0] & 0xFFFF, word[1]);
    }
    printf("[Step 6] Data exchange across columns complete\n");
  }
}

int main(){
  init_virtual_network(Nx, Ny);
  start_virtual_network(&sar_main);
  return 0;
}
