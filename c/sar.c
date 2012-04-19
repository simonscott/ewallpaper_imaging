#include "sar.h"
#include <math.h>
#include <stdio.h>

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
  printf("(%f, %f)", x.real, x.imag);
}

//================================================================================
//==================== Fast Fourier Transform ====================================
//================================================================================

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
void fft_1d_helper(complex* x, int N, int stride, complex* y){
  if(N > 1) {
    fft_1d_helper(x, N/2, 2*stride, y);
    fft_1d_helper(x + stride, N/2, 2*stride, y + N/2);
    
    for(int k=0; k<N/2; k++) {
      complex Wkn1 = c_jexp(-2*pi/N * k);
      complex Wkn2 = c_jexp(-2*pi/N * (k + N/2));
      complex yk1 = c_add(y[k], c_mult(Wkn1, y[k + N/2]));
      complex yk2 = c_add(y[k], c_mult(Wkn2, y[k + N/2]));
      y[k] = yk1;
      y[k+N/2] = yk2;
    }
  } else {
    y[0] = x[0];
  }
}

// fft_1d(complex* x, int N, int stride)
// 1. create a temporary complex buffer of length N
// 2. use fft_1d and compute the fft of x and store it into the temporary buffer
// 3. copy the contents of the temporary buffer to x, with the appropriate stride
void fft_1d(complex* x, int N, int stride){
  complex buffer[N];
  fft_1d(x, N, stride, buffer);
  for(int i=0; i<N; i++)
    x[i*stride] = buffer[i];
}

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
  if(N > 1) {
    ifft_1d_helper(x, N/2, 2*stride, y);
    ifft_1d_helper(x + stride, N/2, 2*stride, y + N/2);
    
    for(int k=0; k<N/2; k++) {
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
  complex buffer[N];
  ifft_1d_helper(x, N, stride, buffer);
  for(int i=0; i<N; i++)
    x[i*stride] = c_scalar_div(buffer[i], N);
}

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
// 7.   compute the interpolated value, x_interp = ratio * x_low + (1 - ratio) * x_high
// 8.   store the interpolated value into the buffer
// 9. copy the values from the temporary buffer back into x
void resample_1d(complex* x, int N, int stride, float* n) {
  complex buffer[N];
  for(int i=0; i<N; i++){
    complex x_interp;
    if(n[i] < 0 || ceil(n[i]) >= N){
      x_interp.real = 0;
      x_interp.imag = 0;
    } else {
      complex x_low = x[(int)floor(n[i]) * stride];
      complex x_high = x[(int)ceil(n[i]) * stride];
      float ratio = n[i] - floor(n[i]);
      x_interp = c_add(c_scalar_mult(x_low, ratio), c_scalar_mult(x_high, 1-ratio));
      buffer[i] = x_interp;
    }
  }

  for(int i=0; i<N; i++)
    x[i*stride] = buffer[i];
}
