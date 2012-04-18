#include "sar.h"
#include <math.h>
#include <stdio.h>

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

inline complex c_add(complex x, complex y) {
  complex r;
  r.real = x.real + y.real;
  r.imag = x.imag + y.imag;
  return r;
}

void fft_1d(complex* x, int N, int stride, complex* y){
  if(N > 1) {
    fft_1d(x, N/2, 2*stride, y); //compute FFT{x_even}
    fft_1d(x + stride, N/2, 2*stride, y + N/2); //compute FFT{x_odd}
    
    for(int k=0; k<N/2; k++) {
      complex Wkn1 = {cos(2*pi/N * k), -sin(2*pi/N * k)};
      complex Wkn2 = {cos(2*pi/N * (k + N/2)), -sin(2*pi/N * (k + N/2))};
      complex yk1 = c_add(y[k], c_mult(Wkn1, y[k + N/2]));
      complex yk2 = c_add(y[k], c_mult(Wkn2, y[k + N/2]));
      y[k] = yk1;
      y[k+N/2] = yk2;
    }
  } else {
    y[0] = x[0];
  }
}


void ifft_1d_helper(complex* x, int N, int stride, complex* y){
  if(N > 1) {
    ifft_1d_helper(x, N/2, 2*stride, y); //compute FFT{x_even}
    ifft_1d_helper(x + stride, N/2, 2*stride, y + N/2); //compute FFT{x_odd}
    
    for(int k=0; k<N/2; k++) {
      complex Wkn1 = {cos(2*pi/N * k), sin(2*pi/N * k)};
      complex Wkn2 = {cos(2*pi/N * (k + N/2)), sin(2*pi/N * (k + N/2))};
      complex yk1 = c_add(y[k], c_mult(Wkn1, y[k + N/2]));
      complex yk2 = c_add(y[k], c_mult(Wkn2, y[k + N/2]));

      y[k] = yk1;
      y[k+N/2] = yk2;
    }
  } else {
    y[0] = x[0];
  }
}

void ifft_1d(complex* x, int N, int stride, complex* y){
  ifft_1d_helper(x, N, stride, y);
  for(int n=0; n<N; n++){
    y[n].real /= N;
    y[n].imag /= N;
  }
}
