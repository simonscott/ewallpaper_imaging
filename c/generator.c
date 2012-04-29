#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

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

void main(){
  filename = "../matlab/small_head.3d";
  Mx = 64;
  My = 64;
  Mz = 113;
  char* scene = read_file(filename, Mx, My, Mz);

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
  
  complex* data = (complex*)malloc(Nx * Ny * Nf * sizeof(complex));
  for(int i=0; i<Nx*Ny*Nf; i++){
    data[i].real = 0;
    data[i].imag = 0;
  }

  int fi;
  #pragma omp parallel shared(scene, data, Nx, Ny, Mx, My, Mz) private(fi)
  {
    #pragma omp for
    for(fi=0; fi<Nf; fi++){
      printf("Frequency %d of %d\n", fi+1, Nf);
      generate_data(scene, data + fi*Nx*Ny, fi, Nx, Ny, Mx, My, Mz);
    }
  }

  write_data("head_scene.dat", data);
}
