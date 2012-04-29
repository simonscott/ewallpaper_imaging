#include <stdlib.h>
#include <stdio.h>
#include "virtual_cpu.h"

// Network Parameters
int nx_cpu;
int ny_cpu;
cpu* cpus;

// Directions
const int up = 0;
const int down = 1;
const int left = 2;
const int right = 3;

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
      args->x = x;
      args->y = y;
      // start the thread
      start_cpu(&(cpus[x + y*nx_cpu]), &forwarding_thread, args);
    }
}

void end_network(){
  for(int i=0; i<nx_cpu*ny_cpu; i++)
    pthread_join(cpus[i].thread, NULL);
}

//================================================================================
//============================ Network Functions =================================
//================================================================================

void send_message(int x, int y, int direction, void* message, int message_size){
  cpu* proc = &cpus[x + y * nx_cpu];

  // up
  if(direction == up){    
    cpu* dest_proc = &cpus[x + (y+1) * nx_cpu];
    send_port(&(dest_proc->down_port), message, message_size);
  }
  // down
  else if(direction == down){
    cpu* dest_proc = &cpus[x + (y-1) * nx_cpu];
    send_port(&(dest_proc->up_port), message, message_size);
  }
  // left
  else if(direction == left){
    cpu* dest_proc = &cpus[(x-1) + y * nx_cpu];
    send_port(&(dest_proc->right_port), message, message_size);
  }
  // right
  else if(direction == right){
    cpu* dest_proc = &cpus[(x+1) + y * nx_cpu];
    send_port(&(dest_proc->left_port), message, message_size);
  }
  // error
  else{
    printf("Unrecognized Direction %d.\n", direction);
    exit(-1);
  }
}

char* receive_message(int x, int y, int direction){
  cpu* proc = &cpus[x + y * nx_cpu];

  // up
  if(direction == up){
    return receive_port(&(proc->up_port));
  }
  // down
  else if(direction == down){
    return receive_port(&(proc->down_port));
  }
  // left
  else if(direction == left){
    return receive_port(&(proc->left_port));
  }
  // right
  else if(direction == right){
    return receive_port(&(proc->right_port));
  }
  // error
  else{
    printf("Unrecognized Direction %d.\n", direction);
    exit(-1);
  }
}

void free_network_port(int x, int y, int direction, void* buffer){
  cpu* proc = &cpus[x + y * nx_cpu];

  // up
  if(direction == up){
    free_port(&(proc->up_port), buffer);
  }
  // down
  else if(direction == down){
    free_port(&(proc->down_port), buffer);
  }
  // left
  else if(direction == left){
    free_port(&(proc->left_port), buffer);
  }
  // right
  else if(direction == right){
    free_port(&(proc->right_port), buffer);
  }
  // error
  else{
    printf("Unrecognized Direction %d.\n", direction);
    exit(-1);
  }  
}

//================================================================================
//============================ Simple Simulation =================================
//================================================================================

void simple_simulation(int x, int y){
  printf("Processor (%d,%d)\n", x, y);
  //Initialize Buffers
  char* up_buffer = (char*)malloc(10*sizeof(int));
  char* down_buffer = (char*)malloc(10*sizeof(int));
  char* left_buffer = (char*)malloc(10*sizeof(int));
  char* right_buffer = (char*)malloc(10*sizeof(int));
  free_network_port(x, y, up, up_buffer);
  free_network_port(x, y, down, down_buffer);
  free_network_port(x, y, left, left_buffer);
  free_network_port(x, y, right, right_buffer);
  //Send message left
  int* msg = (int*)malloc(2*sizeof(int));
  msg[0] = x;
  msg[1] = y;
  if(x > 0)
    send_message(x, y, left, msg, 2*sizeof(int));
  //Receive message right
  if(x < nx_cpu - 1){
    int* r_msg = (int*)receive_message(x, y, right);
    printf("Processor (%d,%d) received [%d %d]\n", x, y, r_msg[0], r_msg[1]);
  }
}

//================================================================================
//=============================== Main Driver ====================================
//================================================================================

void main(){
  init_network(10,10);
  start_network(&simple_simulation);
  end_network();
}
