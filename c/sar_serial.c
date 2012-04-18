#include "sar.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//Initialize the simulation
void init_simulation(complex* data, char* filename){
  FILE* file = fopen(filename, "r");
  
  int n = 0;
  complex num;
  for(int i=0; i < Nx; i++)
    for(int j=0; j < Ny; j++)
      for(int k=0; k < Nf; k++){
        fscanf(file, "%f %f\n", &num.real, &num.imag);
        data[n] = num;
        n++;
      }

  fclose(file);
}

void* safe_malloc(int size, char* error_message){
  void* a = malloc(size);
  if(!a){
    printf("%s", error_message);
    printf("\n");
    exit(-1);
  }
  return a;
}

int main(int argc, char** argv) {

  //Create some fake data
  int N = 16;
  complex x[N];
  for(int i=0; i<N; i++){
    x[i].real = i;
    x[i].imag = 0;
  }

  //Take the FFT
  complex y[N];
  fft_1d(x, N, 1, y);

  //Print out the FFT
  printf("FFT:\n");
  for(int i=0; i<N; i++)
    printf("(%f, %f)\n", y[i].real, y[i].imag);
  printf("\n");

  //Take the IFFT
  complex r[N];
  ifft_1d(y, N, 1, r);

  //Print out the IFFT
  printf("IFFT:\n");
  for(int i=0; i<N; i++)
    printf("(%f, %f)\n", r[i].real, r[i].imag);
  printf("\n");
  
    
  /*
  //Read in data
  complex* s = (complex*)safe_malloc(Nx * Ny * Nf * sizeof(complex),
                         "Failed to allocate memory for radar data.");  
  init_simulation(s, "scene_1.dat");

  //Create temporary buffer space for performing FFTs
  

  //Perform a single 2D FFT for each frequency
  for(int fi = 0; fi < Nf; fi++){
    //Perform 1D FFTs along rows
    for(int i = 0; i < Ny; i++)
      fft_1d(
    //Perform 1D FFTs along columns
  }
  */
}
