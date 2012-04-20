#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include "virtual_network.h"

// Example of two processors bouncing a message around
void* processor_main(int MYTHREAD){
  int N = 100;
  
  if(MYTHREAD == 0){
    int buffer[] = {0};
    send_message(1, (char*)buffer, sizeof(int));
    
    while(1) {
      // Receive message
      int* message = (int*)receive_message(MYTHREAD);
      buffer[0] = *message;
      free_message(MYTHREAD, (char*)message);
      // Process Message
      printf("Thread 0 increments counter from %d to %d\n", buffer[0], buffer[0]+1);
      // Send to Destination
      if(buffer[0] < N){
        buffer[0] = buffer[0] + 1;
        send_message(1, (char*)buffer, sizeof(int));
      }
      else {
        break;
      }
    }
  }
  
  if(MYTHREAD == 1){
    int buffer[1];
    
    while(1){
      // Receive message
      int* message = (int*)receive_message(MYTHREAD);
      buffer[0] = *message;
      free_message(MYTHREAD, (char*)message);
      // Process message
      printf("Thread 1 increments counter from %d to %d\n", buffer[0], buffer[0]+1);
      // Send to Destination
      buffer[0] += 1;
      send_message(0, (char*)buffer, sizeof(int));
      if(buffer[0] >= N)
        break;
    }
  }

  return 0;
}

int main(){
  start_virtual_network(1,2,&processor_main);
  return 0;
}
