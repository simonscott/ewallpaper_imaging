#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>

#define num_processors 10
#define max_messages 100
#define message_memory 1024

//================================================================================
//================== Processors ==================================================
//================================================================================

typedef struct {
  pthread_t thread;
  pthread_cond_t execution;
  pthread_mutex_t execution_lock;

  int num_messages;
  void** messages;
  void* inbox_top;
} processor;

//================================================================================
//====================== Processor Functions =====================================
//================================================================================

processor processors[num_processors];

void* processor_main(void* args);
void init_processor(int i){
  //Synchronization Primitives
  pthread_cond_init(&processors[i].execution, NULL);
  pthread_mutex_init(&processors[i].execution_lock, NULL);

  //Inbox
  processors[i].num_messages = 0;
  processors[i].messages = (void**)malloc(max_messages * sizeof(void*));
  processors[i].inbox_top = malloc(message_memory);

  //Start Thread
  pthread_create(&processors[i].thread, NULL, &processor_main, (void*)i);
}

//Always called with wait_for_messages(MYTHREAD)
void wait_for_messages(int i){
  processor* p = &processors[i];
  pthread_mutex_lock(&p->execution_lock);
  while(p->num_messages == 0)
    pthread_cond_wait(&p->execution, &p->execution_lock);
  pthread_mutex_unlock(&p->execution_lock);
}

void send_message(int i, void* message, int size){
  processor* p = &processors[i];
  //Acquire lock
  pthread_mutex_lock(&p->execution_lock);

  //Add message to inbox
  p->messages[p->num_messages] = p->inbox_top;
  p->num_messages++;
  memcpy(p->inbox_top, message, size);
  p->inbox_top += size;

  //Notify and release lock
  pthread_cond_signal(&p->execution);
  pthread_mutex_unlock(&p->execution_lock);
}

//Always called with receive_message(MYTHREAD)
void* receive_message(int i){
  wait_for_messages(i);
  return processors[i].messages[0];
}

void wait_for_processor(int i){
  pthread_join(processors[i].thread, NULL);
  if(processors[i].num_messages == 0)
    free(processors[i].inbox_top);
  else
    free(processors[i].messages[0]);
}

void free_message(int i){
  processor* p = &processors[i];
  if(p->num_messages == 0){
    printf("Processor %d has an empty inbox.\n", i);
    exit(-1);
  }
  else if(p->num_messages == 1){
    p->inbox_top = p->messages[0];
    p->num_messages = 0;
  }
  else {
    int message_1_size = p->messages[1] - p->messages[0];
    memcpy(p->messages[0], p->messages[1], p->inbox_top - p->messages[1]);
    for(int i=0; i<p->num_messages-1; i++)
      p->messages[i] = p->messages[i+1] - message_1_size;
    p->inbox_top -= message_1_size;
  }
}

//================================================================================
//================================================================================
//================================================================================
int main(){
  printf("Main Driver\n");
  init_processor(0);
  init_processor(1);

  wait_for_processor(0);
  wait_for_processor(1);
  printf("Main Driver Finished\n");

  return 0;
}

void* processor_main(void* args){
  int MYTHREAD = (int)args;
  printf("MYTHREAD = %d\n", MYTHREAD);
  
  if(MYTHREAD == 0){
    printf("Thread 0 waits for thread 1\n");

    int message_out = 42;
    send_message(1, &message_out, sizeof(int));

    int* message = receive_message(MYTHREAD);
    printf("message = %d\n", *message);
    printf("Thread 0 Finished.\n");
  }

  if(MYTHREAD == 1){
    for(int i=0; i<10; i++)
      printf("Thread 1 working: %d\n", i);
    int* message = (int*)receive_message(MYTHREAD);
    printf("THREAD 1: message = %d\n", *message);

    int message_out = -42;
    send_message(0, &message_out, sizeof(int));
    printf("Thread 1 FInished\n");
  }
}
