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
        fscanf(file, "%f, %f\n", &num.real, &num.imag);
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
  // Read in data
  // 1. allocate buffer space for the data. Holds Nx * Ny * Nf complex numbers.
  // 2. pass the filename and buffer to init_simulation which will read the file
  //    into the buffer.
  complex* s = (complex*)safe_malloc(Nx * Ny * Nf * sizeof(complex),
                         "Failed to allocate memory for radar data.");  
  init_simulation(s, "scene_1.dat");  

  // Perform a single 2D FFT for each frequency
  // 1. First do a 1D FFT across each of the Ny rows
  // 2. Then do a 1D FFT across each of the Nx columns
  for(int n = 0; n < Nf; n++){    
    // Perform 1D FFTs along rows
    // 1. Raster order for s is Nx x Ny x Nf
    //    so s[i,j,n] = s[i * Ny * Nf + j * Nf + n]
    // 2. each row starts at (s + j * Nf + n), has length Nx, and stride (Ny * Nf).
    for(int j=0; j<Ny; j++)
      fft_1d_in_place(s + j*Nf + n, Nx, Ny*Nf);

    // Perform 1D FFTs along columns
    // 1. s[i,j,n] = s[i * Ny * Nf + j * Nf + n]
    // 2. each column starts at (s + i*Ny*Nf + n) has length Ny, and stride (Nf).
    for(int i=0; i<Nx; i++)
      fft_1d_in_place(s + i*Ny*Nf + n, Ny, Nf);
  }

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
  //        k  = 2*w/c
  //
  // 3.   compute kz
  //        kz = sqrt( k^2 - kx^2 - ky^2 )
  // 4.   compute the phase delay
  //        phi = exp(j * kz * z0)
  // 5.   multiply the signal with the phase delay
  //        where s(i,j,k) = s[i * Ny * Nf + j * Nf + n]
  for(int i=0; i<Nx; i++)
    for(int j=0; j<Ny; j++)
      for(int n=0; n<Nf; n++){
        float kx = i < Nx/2 ?
          2*pi/Dx * i/Nx :
          2*pi/Dx * (i - Nx)/Nx;
        float ky = j < Ny/2 ?
          2*pi/Dy * j/Ny :
          2*pi/Dy * (j - Ny)/Ny;
        
        float w = 2*pi*(f0 + n*Df);
        float k = 2*w/c_speed;
        float kz = sqrt(k*k - kx*kx - ky*ky);
        
        complex phi = c_jexp(kz * z0);
        s[i * Ny * Nf + j * Nf + n] = c_mult(s[i * Ny * Nf + j * Nf + n], phi);
      }

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
  float w_max = 2*pi * (f0 + (N - 1)*Df);
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
  // 3.   for each step in frequency, n in 0 ... Nf
  //         compute kz = kz_min + (kz_max - kz_min) * n/Nf
  // 4.      compute desired k = 0.5 * sqrt(kx^2 + ky^2 + kz^2)
  // 5.      which corresponds to the interpolated array element
  //            n_interp = (c*k/(2*pi) - f0)/Df
  // 6.      if n_interp lies outside of the range 0 ... Nf, then the interpolated signal is 0
  // 7.      otherwise, compute the interpolated signal :
  // 8.         retrieve the lower element at
  //               s_low = s[i,j,floor(n_interp)] = s[i * Ny * Nf + j * Nf + floor(n_interp)]
  //            retrieve the higher element at
  //               s_high = s[i,j,ceil(n_interp)] = s[i * Ny * Nf + j * Nf + ceil(n_interp)]
  //            interpolate between s_low and s_high
  //               ratio = n_interp - floor(n_interp)
  //               s_interp = ratio * s_low + (1 - ratio) * s_high
  // 9.      store the interpolated signal, s_interp, into s[i,j,n]

  // [ERROR] : We cannot update s in place, we need a temporary buffer
  for(int i=0; i<Nx; i++)
    for(int j=0; j<Ny; j++){
      float kx = i < Nx/2 ?
        2*pi/Dx * i/Nx :
        2*pi/Dx * (i - Nx)/Nx;
      float ky = j < Ny/2 ?
        2*pi/Dy * j/Ny :
        2*pi/Dy * (j - Ny)/Ny;

      for(int n=0; n<Nf; n++){
        float kz = kz_min + (kz_max - kz_min) * n/Nf;
        float k = 0.5 * sqrt(kx*kx + ky*ky + kz*kz);
        float n_interp = (c*k/(2*pi) - f0)/Df;
        
        complex s_interp;
        if(n_interp < 0 || ceil(n_interp) >= Nf) {
          s_interp.real = 0;
          s_interp.imag = 0;
        } else {
          complex s_low = s[i * Ny * Nf + j * Nf + (int)floor(n_interp)];
          complex s_high = s[i * Ny * Nf + j * Nf + (int)ceil(n_interp)];
          float ratio = n_interp - floor(n_interp);
          s_interp.real = ratio * s_low.real + (1 - ratio) * s_high.real;
          s_interp.imag = ratio * s_low.imag + (1 - ratio) * s_high.imag;
          s[i * Ny * Nf + j * Nf + n] = s_interp;
        }
      }
    }

  // Perform an inverse 3D IFFT on the signal
  // 1. First do a 1D IFFT across each of the Ny rows
  // 2. Then do a 1D IFFT across each of the Nx columns
  // 3. Then do a 1D IFFT across each of the Nz lines

  // For each frequency n
  for(int n = 0; n < Nf; n++){    
    // Perform 1D IFFTs along rows
    // 1. Raster order for s is Nx x Ny x Nf
    //    so s[i,j,n] = s[i * Ny * Nf + j * Nf + n]
    // 2. each row starts at (s + j * Nf + n), has length Nx, and stride (Ny * Nf).
    for(int j=0; j<Ny; j++)
      ifft_1d_in_place(s + j*Nf + n, Nx, Ny*Nf);

    // Perform 1D IFFTs along columns
    // 1. s[i,j,n] = s[i * Ny * Nf + j * Nf + n]
    // 2. each column starts at (s + i*Ny*Nf + n) has length Ny, and stride (Nf).
    for(int i=0; i<Nx; i++)
      ifft_1d_in_place(s + i*Ny*Nf + n, Ny, Nf);
  }

  // For each row, i, and column, j
  // 1. s[i,j,n] = s[i * Ny * Nf + j * Nf + n]
  // 2. each line starts at (s + i * Ny * Nf + j * Nf) has length Nf, and stride 1.
  for(int i=0; i<Nx; i++)
    for(int j=0; j<Ny; j++)
      ifft_1d_in_place(s + i * Ny * Nf + j * Nf, Nf, 1);


  // End the simulation by writing out the computed signal and write it out to a file.
  // Pass the computed matrix and a output filename to end_simulation()
  end_simulation(s, "scene_4.out");

  // Free all the temporary memory
  free(s);
}
