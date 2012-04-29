#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <string.h>
#include "virtual_cpu.h"

void init_cpu_port(cpu_port* p){
  // Initialize Message Buffer & Status
  p->buffer = 0;
  p->space_available = 0;
  p->message_waiting = 0;
  // Initialize Notifier
  pthread_cond_init(&(p->signal), NULL);
  pthread_mutex_init(&(p->signal_lock), NULL);
}

void notify_port_listeners(cpu_port* port){
  pthread_mutex_lock(&(port->signal_lock));
  pthread_cond_broadcast(&(port->signal));
  pthread_mutex_unlock(&(port->signal_lock));
}

void send_port(cpu_port* port, void* message, int message_size){
  // Wait for space to become available
  while(!port->space_available){
    pthread_mutex_lock(&(port->signal_lock));
    pthread_cond_wait(&(port->signal), &(port->signal_lock));
    pthread_mutex_unlock(&(port->signal_lock));
  }
  // Copy message to buffer, and notify listeners
  memcpy(port->buffer, (char*)message, message_size);
  port->message_waiting = 1;
  port->space_available = 0;
  notify_port_listeners(port);
}

char* receive_port(cpu_port* port){
  // Wait for message to arrive
  while(!port->message_waiting){
    pthread_mutex_lock(&(port->signal_lock));
    pthread_cond_wait(&(port->signal), &(port->signal_lock));
    pthread_mutex_unlock(&(port->signal_lock));
  }
  // Return message
  return port->buffer;
}

void free_port(cpu_port* port, void* buffer){
  // Set new buffer, and reset state
  port->buffer = (char*)buffer;
  port->message_waiting = 0;
  port->space_available = 1;
  // Notify listeners
  notify_port_listeners(port);
}

void init_cpu(cpu* processor){
  // Initialize Ports
  init_cpu_port(&(processor->left_port));
  init_cpu_port(&(processor->right_port));
  init_cpu_port(&(processor->up_port));
  init_cpu_port(&(processor->down_port));
}

void start_cpu(cpu* processor, void* (*cpu_thread)(void*), void* args){
  int failure = pthread_create(&(processor->thread), NULL, cpu_thread, args);
  if(failure){
    printf("Failed to create processor.\n");
    printf("OS Thread Overload.\n");
    exit(-1);
  }
}

cpu_port* get_port(cpu* processor, int direction){
  if(direction == up)
    return &(processor->up_port);
  else if(direction == down)
    return &(processor->down_port);
  else if(direction == left)
    return &(processor->left_port);
  else if(direction == right)
    return &(processor->right_port);
  else{
    printf("Unrecognized direction: %d\n", direction);
    exit(-1);
  }
}

int opposite_dir(int dir){
  if(dir == up)
    return down;
  else if(dir == down)
    return up;
  else if(dir == left)
    return right;
  else if(dir == right)
    return left;
  else{
    printf("Unrecognized Direction: %d\n", dir);
    exit(-1);
  }
}

void step_in_dir(int dir, int* x, int* y){
  if(dir == up)
    y[0]++;
  else if(dir == down)
    y[0]--;
  else if(dir == left)
    x[0]--;
  else if(dir == right)
    x[0]++;
  else {
    printf("Unrecognized Direction: %d\n", dir);
    exit(-1);
  }
}
