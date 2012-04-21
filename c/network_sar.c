#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include "virtual_network.h"
#include "sar.h"

void dump_data_buffer(char* data_buffer){
  int* id = (int*)data_buffer;
  complex* s = (complex*)(data_buffer + sizeof(int));

  printf("id = %d\n", *id);
  for(int i=0; i<Nf; i++){
    unsigned int* words = (unsigned int*)&s[i];
    printf("[%d, %d, %d]\n", words[0] >> 16, words[0] & 0xFFFF, words[1]);
  }
}

void* sar_main(int MYTHREAD){
  // Get Position
  int xi = MYTHREAD / Ny;
  int yi = MYTHREAD % Ny;

  // Allocate data buffer
  char* data_buffer = (char*)malloc(sizeof(int) + Nf*sizeof(complex));
  complex* s = (complex*)(data_buffer + sizeof(int));
  
  // Generate some fake data
  int* id = (int*)data_buffer;
  id[0] = MYTHREAD;
  for(int i=0; i<Nf; i++){
    int* words = (int*)&s[i];
    words[0] = (xi << 16) + yi;
    words[1] = i;
  }

  if(MYTHREAD == 3)
    dump_data_buffer(data_buffer);

  // Start doing data migration ...
  printf("xi = %d, yi = %d\n", xi, yi);
}

int main(){
  start_virtual_network(Nx, Ny, &sar_main);
  return 0;
}
