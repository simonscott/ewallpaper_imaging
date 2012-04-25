#include "sar.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

//================================================================================
//================ Operations on Complex Numbers =================================
//================================================================================

// Complex Multiplication
// if x is (a,b), and y is (c,d)
// then returns (a*c - b*d, b*c + a*d)
inline complex c_mult(complex x, complex y) {
  complex r;
  float a = x.real;
  float b = x.imag;
  float c = y.real;
  float d = y.imag;
  r.real = a*c - b*d;
  r.imag = b*c + a*d;
  return r;
}

// Complex Addition
inline complex c_add(complex x, complex y) {
  complex r;
  r.real = x.real + y.real;
  r.imag = x.imag + y.imag;
  return r;
}

// Complex Exponential
// Returns exp(j * theta) where theta is specified in radians
inline complex c_jexp(float theta) {
  complex r;
  r.real = cos(theta);
  r.imag = sin(theta);
  return r;
}

// Complex Scalar Multiplication
// Returns (a*x.real, a*x.imag)
inline complex c_scalar_mult(complex x, float a) {
  complex r;
  r.real = x.real * a;
  r.imag = x.imag * a;
  return r;
}

// Complex Scalar Division
// Returns (x.real / a, x.imag / a)
complex c_scalar_div(complex x, float a) {
  complex r;
  r.real = x.real / a;
  r.imag = x.imag / a;
  return r;
}

// c_print(x) :
// Prints out a complex number in form (real, imag)
void c_print(complex x) {
  printf("(%.10f, %.10f)", x.real, x.imag);
}

//================================================================================
//==================== Fast Fourier Transform ====================================
//================================================================================

// Precompute Wkn FFT and IFFT coefficients
// FFT coefficients: Wkn[k] = exp(-2*pi/N * k)
complex* precompute_fft_coefficients(){
  int N = maximum_fft_size;
  complex* Wkn = (complex*)malloc(N * sizeof(complex));
  for(int k=0; k<N; k++)
    Wkn[k] = c_jexp(-2*pi/N * k);
  return Wkn;
}

// IFFT coefficients: Wkn[k] = exp(2*pi/N * k)
// Note that the 1/N scaling factor is not included here.
complex* precompute_ifft_coefficients(){
  int N = maximum_fft_size;
  complex* Wkn = (complex*)malloc(N * sizeof(complex));
  for(int k=0; k<N; k++)
    Wkn[k] = c_jexp(2*pi/N * k);
  return Wkn;
}

// 1D Fourier Transform: fft_1d_helper(x, N, stride, y)
// if N > 1 then
//    1. take the FFT of the even elements in x, and store it in the first half of y
//         x' = x
//         stride' = 2 * stride
//         y' = y
//    2. take the FFT of the odd elements in x, and store it in the second half of y
//         x' = x + stride
//         stride' = 2 * stride
//         y' = y + N/2
//    3. for each frequency step k in 0 ... N/2
//       i. compute the component at Wkn = exp(-j*2*pi/N * k)
//          y[k] = y[k] + Wkn * y[k + N/2]
//       ii. compute the component at Wkn = exp(-j*2*pi/N * (k + N/2))
//          y[k + N/2] = y[k] + Wkn * y[k + N/2]
// otherwise just return the element at x[0]
void fft_1d_helper(complex* x, int N, int stride, complex* y,
                   complex* Wkn, int Wkn_stride){
  if(N > 1) {
    fft_1d_helper(x, N/2, 2*stride, y, Wkn, 2*Wkn_stride);
    fft_1d_helper(x + stride, N/2, 2*stride, y + N/2, Wkn, 2*Wkn_stride);
    
    for(int k=0; k<N/2; k++) {
      complex Wkn1 = Wkn[k*Wkn_stride];
      complex Wkn2 = Wkn[(k + N/2)*Wkn_stride];
      complex yk1 = c_add(y[k], c_mult(Wkn1, y[k + N/2]));
      complex yk2 = c_add(y[k], c_mult(Wkn2, y[k + N/2]));
      y[k] = yk1;
      y[k+N/2] = yk2;
    }
  } else {
    y[0] = x[0];
  }
}

// fft_1d(complex* x, int N, int stride, complex* Wkn)
// 1. create a temporary complex buffer of length N
// 2. use fft_1d and compute the fft of x and store it into the temporary buffer
// 3. copy the contents of the temporary buffer to x, with the appropriate stride
void fft_1d(complex* x, int N, int stride, complex* Wkn){
  // create buffer
  complex buffer[N];

  // fft
  int Wkn_stride = maximum_fft_size / N;
  fft_1d_helper(x, N, stride, buffer, Wkn, Wkn_stride);

  // copy back
  for(int i=0; i<N; i++)
    x[i*stride] = buffer[i];
}

// ifft_1d(complex* x, int N, int stride, complex* Wkn)
// 1. create a temporary complex buffer of length N
// 2. use fft_1d_helper and compute the ifft of x and store it into the temporary buffer
// 3. copy the contents of the temporary buffer to x, with the appropriate stride and
//    divide every element by N (because fft_helper does not do that step).
void ifft_1d(complex* x, int N, int stride, complex* Wkn){
  // create buffer
  complex buffer[N];

  // ifft
  int Wkn_stride = maximum_fft_size / N;  
  fft_1d_helper(x, N, stride, buffer, Wkn, Wkn_stride);

  // copy back
  for(int i=0; i<N; i++)
    x[i*stride] = c_scalar_div(buffer[i], N);
}

/*
//================================================================================
//======================== Inverse Fourier Transform =============================
//================================================================================

// 1D Inverse Fourier Transform: ifft_1d_helper(x, N, stride, y)
// if N > 1 then
//    1. take the IFFT of the even elements in x, and store it in the first half of y
//         x' = x
//         stride' = 2 * stride
//         y' = y
//    2. take the IFFT of the odd elements in x, and store it in the second half of y
//         x' = x + stride
//         stride' = 2 * stride
//         y' = y + N/2
//    3. for each frequency step k in 0 ... N/2
//       i. compute the component at Wkn = exp(j*2*pi/N * k)
//          y[k] = y[k] + Wkn * y[k + N/2]
//       ii. compute the component at Wkn = exp(j*2*pi/N * (k + N/2))
//          y[k + N/2] = y[k] + Wkn * y[k + N/2]
// otherwise just return the element at x[0]
// Note that ifft_1d_helper does *not* do the appropriate 1/N scaling. 
void ifft_1d_helper(complex* x, int N, int stride, complex* y){
  int k;

  if(N > 1) {
    ifft_1d_helper(x, N/2, 2*stride, y);
    ifft_1d_helper(x + stride, N/2, 2*stride, y + N/2);
    
    for(k=0; k<N/2; k++) {
      complex Wkn1 = c_jexp(2*pi/N * k);
      complex Wkn2 = c_jexp(2*pi/N * (k + N/2));
      complex yk1 = c_add(y[k], c_mult(Wkn1, y[k + N/2]));
      complex yk2 = c_add(y[k], c_mult(Wkn2, y[k + N/2]));
      y[k] = yk1;
      y[k+N/2] = yk2;
    }
  } else {
    y[0] = x[0];
  }
}

// ifft_1d(complex* x, int N, int stride)
// 1. create a temporary complex buffer of length N
// 2. use ifft_1d_helper and compute the fft of x and store it into the temporary buffer
// 3. copy the contents of the temporary buffer to x, with the appropriate stride and
//    divide every element by N (because ifft_helper does not do that step).
void ifft_1d(complex* x, int N, int stride){
  int i;
  complex buffer[N];

  ifft_1d_helper(x, N, stride, buffer);

  for(i=0; i<N; i++)
    x[i*stride] = c_scalar_div(buffer[i], N);
}

*/

//================================================================================
//======================= Interpolation ==========================================
//================================================================================

// resample_1d(complex* x, int N, int stride, float* n)
// 1. create a temporary buffer of size N
// 2. for step i in 0 ... N
// 3.   if n[i] is outside of range 0 ... N, then x_interp = 0
// 4.   get the lower element x_low = x[floor(n[i]) * stride]
// 5.   get the higher element x_high = x[ceil(n[i]) * stride]
// 6.   compute the interpolation mixing ratio, ratio = n[i] - floor(n[i])
// 7.   compute the interpolated value, x_interp = (1 - ratio) * x_low + ratio * x_high
// 8.   store the interpolated value into the buffer
// 9. copy the values from the temporary buffer back into x
void resample_1d(complex* x, int N, int stride, float* n) {  
  int i;
  complex buffer[N];

  for(i=0; i<N; i++){
    complex x_interp;
    if(n[i] < 0 || n[i] > N - 1){
      x_interp.real = 0;
      x_interp.imag = 0;
    } else {
      complex x_low = x[(int)floor(n[i]) * stride];
      complex x_high = x[(int)ceil(n[i]) * stride];
      float ratio = n[i] - floor(n[i]);
      x_interp = c_add(c_scalar_mult(x_low, 1 - ratio), c_scalar_mult(x_high, ratio));
    }
    buffer[i] = x_interp;
  }

  for(i=0; i<N; i++)
    x[i*stride] = buffer[i];
}

//================================================================================
//======================= File I/O ===============================================
//================================================================================

// Reads in the file specified by "filename" and stores the data in "data"
// Assumes that the file is formatted as real, imag on each line like so:
// 10.41241, 124.152
// 124.1241, 123.15124
// where the iteration order is Nx * Ny * Nf (Nf being the fastest varying index)
void read_data(complex* data, char* filename){
  FILE* file = fopen(filename, "r");
  if(!file) {
    printf("File %s could not be found.\n", filename);
    exit(-1);
  }
  
  int n = 0;
  int i, j, k;
  complex num;

  for(i=0; i < Nx; i++)
    for(j=0; j < Ny; j++)
      for(k=0; k < Nf; k++){
        fscanf(file, "%f, %f\n", &num.real, &num.imag);
        data[n] = num;
        n++;
      }

  fclose(file);
}

// Outputs the matrix in data to the file specified by "filename".
// Formats the file in the same format that it was written in as per above.
void write_data(complex* data, char* filename){
  int i, j, k;
  FILE* file = fopen(filename, "w");
  int n = 0;

  for(i=0; i<Nx; i++)
    for(j=0; j<Ny; j++)
      for(k=0; k<Nf; k++){
        fprintf(file, "%f, %f\n", data[n].real, data[n].imag);
        n++;
      }
  
  fclose(file);
}

//================================================================================
//======================= Memory Allocation / Deallocation =======================
//================================================================================

// Safely allocate a block of memory
// Internally calls malloc to allocate memory
// If malloc fails, then error_message is printed out and the program is aborted.
void* safe_malloc(int size, char* error_message) {
  void* a = malloc(size);
  if(!a){
    printf("%s", error_message);
    printf("\n");
    exit(-1);
  }
  return a;
}

//================================================================================
//======================= Timing Functions =======================================
//================================================================================

// on tick() save the current_time to tick_time
struct timeval tick_time;

void tick(){
  gettimeofday(&tick_time, NULL);
  printf("Tick\n");
}

void tock(){
  struct timeval tock_time;
  gettimeofday(&tock_time, NULL);
  long seconds = tock_time.tv_sec - tick_time.tv_sec;
  long milliseconds = tock_time.tv_usec - tick_time.tv_usec;
  double time = seconds + milliseconds/1e6;
  printf("Tock: Time Passed = %f\n", time);
}
