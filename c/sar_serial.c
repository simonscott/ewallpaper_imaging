#include "sar.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

// Precomputed fft and ifft coefficients
complex* Wkn_fft;
complex* Wkn_ifft;

// Perform a single 2D FFT by performing nx+ny 1D FFTs
// fft_2d(complex* x, int nx, int ny, int x_stride, int y_stride)
// 1. Perform 1D FFTs along rows
//    Element x[i,j] is located at x[i * x_stride + j * y_stride]
//    So each row starts at x + 0 * x_stride + j * y_stride, has length nx, and stride x_stride
// 2. Perform 1D FFTs along columns
//    Element x[i,j] is located at x[i * x_stride + j * y_stride]
//    So each column starts at x + i * x_stride + 0 * y_stride, has length ny, and stride y_stride
void fft_2d(complex* x, int nx, int ny, int x_stride, int y_stride) {
  int i, j;

  for(j=0; j<ny; j++)
    fft_1d(x + j * y_stride, nx, x_stride, Wkn_fft);
  for(i=0; i<nx; i++)
    fft_1d(x + i * x_stride, ny, y_stride, Wkn_fft);
}

// Perform a single 2D IFFT by performing nx+ny 1D IFFTs
// ifft_2d(complex* x, int nx, int ny, int x_stride, int y_stride)
// 1. Perform 1D IFFTs along rows
//    Element x[i,j] is located at x[i * x_stride + j * y_stride]
//    So each row starts at x + 0 * x_stride + j * y_stride, has length nx, and stride x_stride
// 2. Perform 1D IFFTs along columns
//    Element x[i,j] is located at x[i * x_stride + j * y_stride]
//    So each column starts at x + i * x_stride + 0 * y_stride, has length ny, and stride y_stride
void ifft_2d(complex* x, int nx, int ny, int x_stride, int y_stride) {
  int i, j;

  for(j=0; j<ny; j++)
    ifft_1d(x + j * y_stride, nx, x_stride, Wkn_ifft);
  for(i=0; i<nx; i++)
    ifft_1d(x + i * x_stride, ny, y_stride, Wkn_ifft);
}

// Perform a single 3D IFFT by performing NZ 2D IFFTs followed by nx*ny 1D IFFTs
// 1. Perform 2D IFFTs along xy planes
//    Element x[i,j] is located at x[i * x_stride + j * y_stride + k * z_stride]
//    So each plane starts at x + 0 * x_stride + 0 * y_stride + k * z_stride
// 2. Perform 1D IFFTs along z axis
//    Each line starts at x + i * x_stride + j * y_stride + 0 * z_stride
void ifft_3d(complex* x, int nx, int ny, int nz, int x_stride, int y_stride, int z_stride) {
  int i, j, k;

  for(k=0; k<nz; k++)
    ifft_2d(x + k * z_stride, nx, ny, x_stride, y_stride);

  for(i=0; i<nx; i++)
    for(j=0; j<ny; j++)
      ifft_1d(x + i * x_stride + j * y_stride, nz, z_stride, Wkn_ifft);
}


void read_original_data(complex* data, char* filename){
  FILE* file = fopen(filename, "r");
  if(!file) {
    printf("File %s could not be found.\n", filename);
    exit(-1);
  }
  
  int n = 0;
  int i, j, k;
  complex num;

  for(i=0; i < 128; i++)
    for(j=0; j < 128; j++)
      for(k=0; k < 256; k++){
        fscanf(file, "%f, %f\n", &num.real, &num.imag);
        if(i % (128/Nx) == 0 &&
           j % (128/Ny) == 0 &&
           k % (256/Nf) == 0){
          data[n] = num;
          n++;
        }
      }

  fclose(file);
}

int main(int argc, char** argv) {
  //Precompute FFT coefficients
  Wkn_fft = precompute_fft_coefficients();
  Wkn_ifft = precompute_ifft_coefficients();
  
  // Declare local variables
  int i, j, n;

  // Read in data
  // 1. allocate buffer space for the data. Holds Nx * Ny * Nf complex numbers.
  // 2. pass the filename and buffer to read_data which will read the file
  //    into the buffer.
  printf("Reading Data ...\n");
  tick();
  complex* s = (complex*)safe_malloc(Nx * Ny * Nf * sizeof(complex),
                         "Failed to allocate memory for radar data.");
  read_original_data(s, "scene_4.dat");
  //read_data(s, "scene_4.dat");
  tock();

  // Perform a single 2D FFT for each frequency
  // Each element s[i,j,n] is located at s[i * Ny * Nf + j * Nf + n]
  // Thus x-stride = Ny*Nf, y-stride = Nf, and z-stride = 1
  // and each xy plane starts at s + 0 * Ny * Nf + 0 * Nf + n
  printf("Performing FFT\n");
  tick();
  for(n=0; n<Nf; n++)
    fft_2d(s + n, Nx, Ny, Ny*Nf, Nf);
  tock();

  // Multiply each element in the frequency-domain signal by the
  // downward continuation phase operator.
  // 1. for each element (i,j,n) in the signal matrix
  //      i in 0 ... Nx, j in 0 ... Ny, n in 0 ... Nf
  // 2.   compute kx, ky, and k
  //        kx = 2*pi/Dx * i/Nx          if i < Nx/2,
  //             2*pi/Dx * (i-Nx)/Nx     otherwise
  //        ky = 2*pi/Dy * j/Ny          if j < Ny/2,
  //             2*pi/Dy * (j-Ny)/Ny     otherwise
  //        w  = 2*pi*(f0 + n*Df)
  //        k  = w/c
  //
  // 3.   compute kz
  //        kz = sqrt(4 * k^2 - kx^2 - ky^2 )
  // 4.   compute the phase delay
  //        phi = exp(j * kz * z0)
  // 5.   multiply the signal with the phase delay
  //        where s(i,j,k) = s[i * Ny * Nf + j * Nf + n]
  printf("Performing Downward Continuation.\n");
  tick();
  for(i=0; i<Nx; i++)
    for(j=0; j<Ny; j++)
      for(n=0; n<Nf; n++){
        float kx = i < Nx/2 ?
          2*pi/Dx * i/Nx :
          2*pi/Dx * (i - Nx)/Nx;
        float ky = j < Ny/2 ?
          2*pi/Dy * j/Ny :
          2*pi/Dy * (j - Ny)/Ny;
        
        float w = 2*pi*(f0 + n*Df);
        float k = w/c_speed;
        float kz = sqrt(4*k*k - kx*kx - ky*ky);
        
        complex phi = c_jexp(kz * z0);
        s[i * Ny * Nf + j * Nf + n] = c_mult(s[i * Ny * Nf + j * Nf + n], phi);        
      }
  tock();

  // Calculate the range of the Stolt interpolation indices.
  // The minimum angular frequency, w_min = 2*pi * f0
  // The maximum angular frequency, w_max = 2*pi * (f0 + (N - 1)*Df)
  // From which the
  //   minimum wavenumber, k_min = w_min / c
  //   maximum wavenumber, k_max = w_max / c
  // The maximum wavenumber in the x direction, kx_max = 2*pi/Dx * 0.5 * (Nx-1)/Nx
  // The maximum wavenumber in the y direction, ky_max = 2*pi/Dy * 0.5 * (Ny-1)/Ny
  // The minimum wavenumbers in the x and y direction are assumed to be 0
  // From which the
  //   minimum wavenumber in the z direction, kz_min = sqrt(4*k_min^2 - kx_max^2 - ky_max^2)
  //   maximum wavenumber in the z direction, kz_max = sqrt(4*k_max^2 - 0 - 0)
  float w_min = 2*pi * f0;
  float w_max = 2*pi * (f0 + (Nf - 1)*Df);
  float k_min = w_min / c_speed;
  float k_max = w_max / c_speed;
  float kx_max = 2*pi/Dx * 0.5 * (Nx-1)/Nx;
  float ky_max = 2*pi/Dy * 0.5 * (Ny-1)/Ny;
  float kz_min = sqrt(4*k_min*k_min - kx_max*kx_max - ky_max*ky_max);
  float kz_max = sqrt(4*k_max*k_max);

  // Perform Stolt Interpolation
  // 1. for each step in the x direction, i in 0 ... Nx
  //    and each step in the y direction, j in 0 ... Ny
  // 2.   compute kx, and ky as per step 2. above
  // 3.   create float buffer of size Nf for storing the interpolation indices, n_interp
  // 4.   for each step in frequency, n in 0 ... Nf
  //         compute kz = kz_min + (kz_max - kz_min) * n/(Nf - 1)
  // 4.      compute desired k = 0.5 * sqrt(kx^2 + ky^2 + kz^2)
  // 5.      which corresponds to the interpolated array element
  //            n_interp[n] = (c*k/(2*pi) - f0)/Df
  // 6.   resample this line in s on interpolated indices n_interp
  //         s[i,j,n] is at s[i * Ny * Nf + j * Nf + n] thus this line
  //           starts at s + i * Ny * Nf + j * Nf + 0, has length Nf, and has stride 1
  printf("Performing Stolt Interpolation.\n");
  tick();
  for(i=0; i<Nx; i++)
    for(j=0; j<Ny; j++){
      float kx = i < Nx/2 ?
        2*pi/Dx * i/Nx :
        2*pi/Dx * (i - Nx)/Nx;
      float ky = j < Ny/2 ?
        2*pi/Dy * j/Ny :
        2*pi/Dy * (j - Ny)/Ny;

      float n_interp[Nf];
      for(n=0; n<Nf; n++){
        float kz = kz_min + (kz_max - kz_min) * n/(Nf - 1);
        float k = 0.5 * sqrt(kx*kx + ky*ky + kz*kz);
        n_interp[n] = (c_speed*k/(2*pi) - f0)/Df;
      }
      
      resample_1d(s + i*Ny*Nf + j*Nf, Nf, 1, n_interp);
    }
  tock();

  // Perform a 3D IFFT on the signal
  // Each element s[i,j,n] is located at s[i * Ny * Nf + j * Nf + n]
  // Thus x-stride = Ny*Nf, y-stride = Nf, and z-stride = 1
  printf("Performing IFFT.\n");
  tick();
  ifft_3d(s, Nx, Ny, Nf, Ny*Nf, Nf, 1);
  tock();

  // End the simulation by writing out the computed signal and write it out to a file.
  // Pass the computed matrix and a output filename to write_data()
  printf("Writing data ...\n");
  tick();
  write_data(s, "scene_4.out");
  tock();
  printf("Done.\n");

  // Free all the temporary memory
  free(s);
}
