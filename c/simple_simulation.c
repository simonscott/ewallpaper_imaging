#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "network.h"
#include "virtual_cpu.h"
#include "sar.h"

//================================================================================
//============================ Simple Simulation =================================
//================================================================================

void get_simulation_size(int* x, int* y){
  x[0] = Nx;
  y[0] = Ny;
}

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

char* receive_line(receive_status* status, int less_dir){
  int more_dir = opposite_dir(less_dir);
  
  if(status->num_less == 0 ||
     (status->num_more > 0 && status->receive_dir == more_dir)){
      status->receive_dir = less_dir;
      // Receive and copy message into buffer
      char* msg = receive_message(status->x, status->y, more_dir);
      // Forward message
      memcpy(status->buffer, (char*)msg, status->size);
      free_network_port(status->x, status->y, more_dir, msg);
      if(status->x > 0)
        send_message(status->x, status->y, less_dir, status->buffer, status->size);
      // Decrement counter
      status->num_more--;
      // Process message
      return status->buffer;
  }
  else if(status->num_more == 0 ||
          (status->num_less > 0 && status->receive_dir == less_dir)){
    status->receive_dir = more_dir;
    // Receive and copy message into buffer
    char* msg = receive_message(status->x, status->y, less_dir);
    // Forward message
    memcpy(status->buffer, (char*)msg, status->size);
    free_network_port(status->x, status->y, less_dir, msg);
    if(status->x < Nx - 1)
      send_message(status->x, status->y, more_dir, status->buffer, status->size);
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
  return receive_line(status, left);
}

char* receive_col(receive_status* status){
  return receive_line(status, down);
}

void network_simulation(int x, int y){
  //Initialize Buffers
  int msg_size = Nf * sizeof(complex) + 2 * sizeof(int);
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
  complex* s = (complex*)(data + 2*sizeof(int));
  for(int i=0; i<Nf; i++){
    int* words = (int*)(&s[i]);
    words[0] = (x << 16) + y;
    words[1] = i;
  }

  //Send row
  send_row(x, y, data, msg_size);

  //Receive row
  char* buffer = (char*)malloc(msg_size);
  receive_status status = start_receive_row(x, y, buffer, msg_size);
  for(int i=0; i<Nx-1; i++){
    int* msg = (int*)receive_line(&status, left);
    printf("Processor (%d,%d) received [%d,%d]\n", x, y, msg[0], msg[1]);
  }
}
