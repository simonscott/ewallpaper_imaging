#include <stdlib.h>
#include <stdio.h>
#include "virtual_cpu.h"
#include "network.h"

// Network Parameters
int nx_cpu;
int ny_cpu;
cpu* cpus;

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
//=============================== Main Driver ====================================
//================================================================================

// Provided Simulation
void get_simulation_size(int* num_x, int* num_y);
void network_simulation(int x, int y);

void main(){
  int num_x, num_y;
  get_simulation_size(&num_x, &num_y);
  init_network(num_x,num_y);
  start_network(&network_simulation);
  end_network();
}
