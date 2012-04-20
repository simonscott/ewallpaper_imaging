#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include "virtual_network.h"

void* processor_main(int MYTHREAD){
  printf("MYTHREAD = %d\n", MYTHREAD);
  
  if(MYTHREAD == 0){
    int message_out = 42;
    send_message(1, &message_out, sizeof(int));
    printf("Thread 0 waits for thread 1\n");
    
    int* message = (int*)receive_message(MYTHREAD);
    printf("message = %d\n", *message);
    
    printf("Thread 0 Finished.\n");
  }

  if(MYTHREAD == 1){
    for(int i=0; i<10; i++)
      printf("Thread 1 working: %d\n", i);
    
    int* message = (int*)receive_message(MYTHREAD);
    printf("THREAD 1: message = %d\n", *message);
    free_message(MYTHREAD, message);

    int message_out = -42;
    send_message(0, &message_out, sizeof(int));
    printf("Thread 1 FInished\n");
  }

  return 0;
}

int main(){
  start_virtual_network(1,2,&processor_main);
  return 0;
}
