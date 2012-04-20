#include <stdlib.h>
#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include "virtual_network.h"

//================================================================================
//==================== Network Parameters ========================================
//================================================================================

const int max_messages = 100;
const int message_memory = 1024;
int rows_of_processors;
int cols_of_processors;
int num_processors;

//================================================================================
//===================== Processors ===============================================
//================================================================================

typedef struct {
  pthread_t thread;
  pthread_cond_t execution;
  pthread_mutex_t execution_lock;

  int num_messages;
  void** messages;
  void* inbox_top;
} processor;

processor* processors;

//================================================================================
//======================== Internal Processor Functions ==========================
//================================================================================

void init_processor(int i, processor_main_function main_function){
  //Synchronization Primitives
  pthread_cond_init(&processors[i].execution, NULL);
  pthread_mutex_init(&processors[i].execution_lock, NULL);

  //Inbox
  processors[i].num_messages = 0;
  processors[i].messages = (void**)malloc(max_messages * sizeof(void*));
  processors[i].inbox_top = malloc(message_memory);

  //Start Thread
  void* (*thread_entry)(void*) = (void* (*)(void*))main_function;
  pthread_create(&processors[i].thread, NULL, thread_entry, (void*)i);
}

void wait_for_messages(processor* p){
  pthread_mutex_lock(&(p->execution_lock));
  while(p->num_messages == 0)
    pthread_cond_wait(&(p->execution), &(p->execution_lock));
  pthread_mutex_unlock(&(p->execution_lock));
}

void add_message(int i, void* message, int size){
  processor* p = &processors[i];
  if(p->num_messages >= max_messages){
    printf("Maximum number of messages exceeded on processor %d\n", i);
    exit(-1);
  }

  if(p->num_messages > 0){
    long end_of_memory = (long)p->messages[0] + message_memory;
    long requested_end = (long)p->inbox_top + size;
    if(requested_end >= end_of_memory){
      printf("Message buffer memory exceeded on processor %d\n", i);
      exit(-1);
    }
  }
  
  p->messages[p->num_messages] = p->inbox_top;
  p->num_messages++;
  memcpy(p->inbox_top, message, size);
  p->inbox_top = (void*)((long)p->inbox_top + size);  
}

void send_message(int i, void* message, int size){
  processor* p = &processors[i];
  //Acquire lock
  pthread_mutex_lock(&(p->execution_lock));

  //Add message to inbox
  add_message(i, message, size);

  //Notify and release lock
  pthread_cond_signal(&(p->execution));
  pthread_mutex_unlock(&(p->execution_lock));
}

void* receive_message(int threadid){
  wait_for_messages(&processors[threadid]);
  return processors[threadid].messages[0];
}

void free_message(int i, void* message){
  processor* p = &processors[i];
  
  if(p->num_messages == 0){
    printf("Processor %d has an empty inbox.\n", i);
    exit(-1);
  }
  
  if(p->messages[0] != message){
    printf("You must free messages in order.\n");
    exit(-1);
  }
  
  if(p->num_messages == 1){
    p->inbox_top = p->messages[0];
    p->num_messages = 0;
  }
  else {
    int copy_size = (long)p->inbox_top - (long)p->messages[1];
    memcpy(p->messages[0], p->messages[1], copy_size);
    
    int message_1_size = (long)p->messages[1] - (long)p->messages[0];
    for(int i=0; i<p->num_messages-1; i++)
      p->messages[i] = (void*)((long)p->messages[i+1] - message_1_size);
    p->inbox_top = (void*)((long)p->inbox_top - message_1_size);
  }
}

void free_processor(processor* p){
  if(p->num_messages == 0)
    free(p->inbox_top);
  else
    free(p->messages[0]);
}

//================================================================================
//===================== Network Functions ========================================
//================================================================================
void start_virtual_network(int rows, int cols, processor_main_function p_main){
  //Set Parameters
  rows_of_processors = rows;
  cols_of_processors = cols;
  num_processors = rows_of_processors * cols_of_processors;

  //Allocate Processors
  processors = (processor*)malloc(num_processors * sizeof(processor));

  //Initialize Processors
  for(int i=0; i<num_processors; i++)
    init_processor(i, p_main);

  //Wait for Processors to finish processing
  for(int i=0; i<num_processors; i++)
    pthread_join(processors[i].thread, NULL);

  //Free Processor Memory
  for(int i=0; i<num_processors; i++)
    free_processor(&processors[i]);

  //Free Processors
  free(processors);
}
