#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <pthread.h>

#define M_PI 3.14159265358979323846

typedef struct {
  float real;
  float imag;
} complex;

char* filename;
int Nx, Ny, Nf;
int Mx, My, Mz;

float f0, Df;
float Dx, Dy, Dz;
float z0;

char* read_file(char* filename, int Mx, int My, int Mz){
  FILE* file = fopen(filename, "r");
  if(!file){
    printf("%s could not be found!\n", filename);
    exit(-1);
  }
  
  int size = Mx * My * Mz;
  char* ret = (char*)malloc(size);
  fread(ret, 1, size, file);
  fclose(file);
  return ret;
}

void generate_data(char* scene, complex* data, int fi,
                   int Nx, int Ny,
                   int Mx, int My, int Mz){  
  for(int i=0; i < Mx; i++){
    for(int j=0; j < My; j++){
      for(int k=0; k < Mz; k++){
        float r = scene[i + j * Mx + k * Mx * My];
        if(r != 0){
          float tx = (i - (Mx - 1.0) / 2.0) * Dx;
          float ty = (j - (My - 1.0) / 2.0) * Dy;
          float tz = z0 + k * Dz;
          for(int aj=0; aj<Ny; aj++){
            for(int ai=0; ai<Nx; ai++){
              float ax = (ai - (Nx - 1.0) / 2.0) * Dx;
              float ay = (aj - (Ny - 1.0) / 2.0) * Dy;
              float f = f0 + fi * Df;
              float d = sqrt((tx-ax)*(tx-ax) + (ty-ay)*(ty-ay) + tz*tz);
              float c_speed = 299792458;
              float dt = 2*d/c_speed;
              float phi = dt * 2 * M_PI * f;
              float a_real = r * cos(phi);
              float a_imag = -r * sin(phi);
              data[ai + aj*Nx].real += a_real;
              data[ai + aj*Nx].imag += a_imag;
            }
          }        
        }
      }
    }
  }
}

void write_data(char* filename, complex* data){
  FILE* file = fopen(filename, "w");

  for(int i=0; i<Nx; i++)
    for(int j=0; j<Ny; j++)
      for(int k=0; k<Nf; k++){
        complex c = data[i + j*Nx + k*Nx*Ny];
        fprintf(file, "%f, %f\n", c.real, c.imag);
      }
  
  fclose(file);
}

// Shared Thread Variables
char* scene;
complex* data;

void* generate_thread(void* args){
  int* f_args = (int*)args;
  int f_min = f_args[0];
  int f_max = f_args[1];
  for(int fi = f_min; fi < f_max; fi++){
    printf("Frequency %d\n", fi);
    generate_data(scene, data + fi*Nx*Ny, fi, Nx, Ny, Mx, My, Mz);
  }
  return 0;
}

void main(){
  filename = "../matlab/small_head.3d";
  Mx = 64;
  My = 64;
  Mz = 113;
  printf("Reading File %s\n", filename);
  scene = read_file(filename, Mx, My, Mz);

  Nx = 128;
  Ny = 128;
  Nf = 256;
  z0 = 1;
  f0 = 10e9;
  Dx = 0.012;
  Dy = 0.012;
  Dz = 0.012;
  float B = 2e9;
  Df = B/Nf * 256/Nf;
  
  data = (complex*)malloc(Nx * Ny * Nf * sizeof(complex));
  for(int i=0; i<Nx*Ny*Nf; i++){
    data[i].real = 0;
    data[i].imag = 0;
  }

  //Launch generator threads
  int n_threads = 24;
  pthread_t threads[n_threads];
  int thread_args[n_threads][2];

  //Compute thread work
  for(int i=0; i<n_threads; i++){
    thread_args[i][0] = i * Nf / n_threads;
    thread_args[i][1] = (i + 1) * Nf / n_threads;
    if(thread_args[i][1] > Nf)
      thread_args[i][1] = Nf;
  }

  //Start threads
  printf("Starting Computation\n");
  for(int i=0; i<n_threads; i++){
    int failure = pthread_create(&threads[i], NULL, &generate_thread, &thread_args[i][0]);
    if(failure){
      printf("Failed to create thread\n");
      exit(-1);
    }
  }

  //End threads
  for(int i=0; i<n_threads; i++)
    pthread_join(threads[i], NULL);

  write_data("small_head.dat", data);
}
