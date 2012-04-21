#include <stdlib.h>
#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include "sar.h"
#include "virtual_network.h"

//================================================================================
//==================== Network Parameters ========================================
//================================================================================

// max_messages holds the maximum number of messages that a single node can hold.
// message_memory holds the amount of memory in the network buffer per node.
// rows_of_processors and cols_of_processors indicate the number of rows and
//   columns of processors in the network.
// num_processors = rows_of_processors * cols_of_processors is calculated upon
//   initialization.

const int max_messages = 100;
//const int message_memory = (Nf + 16)*sizeof(complex);
int rows_of_processors;
int cols_of_processors;
int num_processors;

//================================================================================
//===================== Processors ===============================================
//================================================================================

// Processor struct.
// Contains fields, thread, execution, and execution_lock for managing concurrency.
// Contains fields num_messages, messages, and inbox_top for managing the message inbox.
// num_messages holds the number of messages currently queued in the message inbox.
// messages is a char* array where each element holds a pointer the beginning of the
//    message
// inbox_top is a pointer to the next available byte in the message inbox.
typedef struct {
  pthread_t thread;
  pthread_cond_t execution;
  pthread_mutex_t execution_lock;

  int num_messages;
  char** messages;
  char* inbox_top;
} processor;

// All the processors in the simulation are stored here.
// Thus processor 42 refers to processors[42]
processor* processors;

//================================================================================
//======================== Internal Processor Functions ==========================
//================================================================================

// init_processor(i) :
// 1. Initialize synchronization primitives, execution and execution_lock
// 2. Initialize the message inbox, 0 messages, space to hold pointers to
//    max_messages number of messages, and space to hold messages
void init_processor(int i){
  //Synchronization Primitives
  pthread_cond_init(&processors[i].execution, NULL);
  pthread_mutex_init(&processors[i].execution_lock, NULL);

  //Inbox
  processors[i].num_messages = 0;
  processors[i].messages = (char**)malloc(max_messages * sizeof(char*));
  processors[i].inbox_top = (char*)malloc(message_memory);
}

// start_processor(i, main_function) :
// 1. Creates the thread for processor i, and starts it.
// Note that this function should not be called until all the other processors
// that this processor might talk to are already initialized.
void start_processor(int i, processor_main_function main_function){
  void* (*thread_entry)(void*) = (void* (*)(void*))main_function;
  pthread_create(&processors[i].thread, NULL, thread_entry, (void*)(long)i);
}

// wait_for_messages(p) :
// Blocks until there is a message in the message inbox.
// Called internally by receive_virtual_message()
// 1. Acquire the execution_lock
// 2. While there are no number of messages,
// 3. Wait on the condition variable. (Expects the condition variable to be
//    notified on arrival messages).
// 4. Release the execution_lock
void wait_for_messages(processor* p){
  pthread_mutex_lock(&(p->execution_lock));
  while(p->num_messages == 0)
    pthread_cond_wait(&(p->execution), &(p->execution_lock));
  pthread_mutex_unlock(&(p->execution_lock));
}

// add_message(i, message, size) :
// Adds the specified size to the message inbox of processor i.
// Will be called internally by send_virtual_message.
// 1. Retrieve the appropriate processor, p
// 2. Print an error if the message inbox is full.
// 3. Print an error if the message inbox is out of memory
// Copy the message to the buffer:
//   4. Copy the message to the message buffer at inbox_top
//   5. Store a pointer to the message in the inbox
//   6. Increase the number of messages by 1
//   7. Increase inbox_top by the size of the message
void add_message(int i, char* message, int size){
  //Retrieve processor
  processor* p = &processors[i];
  //Maximum number of messages reached
  if(p->num_messages >= max_messages){
    printf("Maximum number of messages exceeded on processor %d\n", i);
    exit(-1);
  }
  //Buffer overflow
  if(p->num_messages > 0){
    char* end_of_memory = p->messages[0] + message_memory;
    if(p->inbox_top + size >= end_of_memory){
      printf("Message buffer memory exceeded on processor %d\n", i);
      exit(-1);
    }
  }
  //Copy message to buffer
  memcpy(p->inbox_top, message, size);
  p->messages[p->num_messages] = p->inbox_top;
  p->num_messages++;
  p->inbox_top += size;
}

// send_virtual_message(i, message, size)
// Sends a message to the processor i. Copies the given message to the message buffer
// of the destination processor.
// 1. Retrieve the processor, p
// 2. Acquire the execution lock
// 3. Add the message the inbox
// 4. Notify the processor's the condition variable
// 3. Release the execution lock
void send_virtual_message(int i, char* message, int size){
  processor* p = &processors[i];
  pthread_mutex_lock(&(p->execution_lock));
  add_message(i, message, size);
  pthread_cond_signal(&(p->execution));
  pthread_mutex_unlock(&(p->execution_lock));
}

// receive_virtual_message(MYTHREAD) :
// Returns a pointer to the first message in the message inbox
// If the message inbox is empty, then wait for a message to arrive.
char* receive_virtual_message(int threadid){
  if(processors[threadid].num_messages == 0)
    wait_for_messages(&processors[threadid]);
  return processors[threadid].messages[0];
}

// free_virtual_message(MYTHREAD, message) :
// Frees the first message in the inbox (which is pointed to by message).
// 1. Retrieve the processor
// 2. Print an error if the inbox is empty and there is no message to free.
// 3. Print an error if the first message in the inbox is not the same as message.
// 4. If there is only one message in the inbox :
//       then message[0] just holds a pointer to the start of the message memory
// 5.    so just reset inbox_top to point to message[0]
// 6.    and reset num_messages to 0.
// 7. If there is more than one message :
// 8.    Copy messages[1 ... end] to messages[0]
// 9.    Compute the size of the first message, message_1_size
// 10.   Decrease the inbox_top by message_1_size
// 10.   For each message[i], update its pointer to be (what it used to be) - message_1_size
void free_virtual_message(int i, char* message){
  // Retrieve Processor
  processor* p = &processors[i];
  // Inbox Empty
  if(p->num_messages == 0){
    printf("Processor %d has an empty inbox.\n", i);
    exit(-1);
  }
  // Given messages doesn't match message[0]
  if(p->messages[0] != message){
    printf("You must free messages in order.\n");
    exit(-1);
  }
  // One message in inbox
  if(p->num_messages == 1){
    p->inbox_top = p->messages[0];
    p->num_messages = 0;
  }
  // Multiple messages in inbox
  else {
    memcpy(p->messages[0], p->messages[1], p->inbox_top - p->messages[1]);    
    int message_1_size = p->messages[1] - p->messages[0];
    for(int i=0; i < p->num_messages-1; i++)
      p->messages[i] = p->messages[i+1] - message_1_size;
    p->inbox_top -= message_1_size;
  }
}

// free_processor(p) :
// Free the memory allocated for the message inbox.
// If message inbox is empty, then inbox_top points to the start of the inbox
// Otherwise, messages[0] points to the start of the inbox
void free_processor(processor* p){
  char* message_buffer = p->num_messages == 0 ?
    p->inbox_top :
    p->messages[0];
  free(message_buffer);
}

//================================================================================
//===================== Network Functions ========================================
//================================================================================
// start_virtual_network(rows, cols, processor_main)
// Starts the virtual network with the given configuration, runs the simulation,
// and performs cleanup and returns.
// 1. Set the global parameters of the network
// 2. Allocate the processors
// 3. Initialize all the processors
// 4. Start all the processors
// 5. Wait for all the processor threads to finish
// 6. Free all the processors, and the global processors array
void start_virtual_network(int rows, int cols, processor_main_function p_main){
  //Set Parameters
  rows_of_processors = rows;
  cols_of_processors = cols;
  num_processors = rows_of_processors * cols_of_processors;
  //Allocate Processors
  processors = (processor*)malloc(num_processors * sizeof(processor));
  //Initialize Processors
  for(int i=0; i<num_processors; i++)
    init_processor(i);
  //Start Processors
  for(int i=0; i<num_processors; i++)
    start_processor(i, p_main);
  //Wait for Processors to finish processing
  for(int i=0; i<num_processors; i++)
    pthread_join(processors[i].thread, NULL);
  //Free Processor Memory and Processor Array
  for(int i=0; i<num_processors; i++)
    free_processor(&processors[i]);
  free(processors);
}
