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

void fft_helper(complex* x, int N, int stride, complex* y){
  if(N > 1) {
    
  } else {
    y[0] = x[0];
  }
}

void fft_1d(complex* x, int N, int stride) {

  if(N > 1){
    fft_1d(x, N/2, stride*2); //compute FFT{x_even}
    
    fft_1d(x+stride, N/2, stride*2); //compute FFT{x_odd}


    printf("Even FFT: ");
    for(int i=0; i<N/2; i++)
      printf("(%f,%f), ", x[i*stride*2].real, x[i*stride*2].imag);
    printf("\n");

    printf("Odd FFT: ");
    for(int i=0; i<N/2; i++)
      printf("(%f,%f), ", x[i*stride*2 + stride].real, x[i*stride*2 + stride].imag);
    printf("\n");
    
    for(int k=0; k<N/2; k++){
      //Positive Freq
      complex Xk1;
      complex Xk2;
      
      complex Wkn1 = {cos(2*pi/N * k), -sin(2*pi/N * k)};
      printf("(%f,%f) + (%f,%f)*(%f,%f)\n",
             x[2*k*stride].real, x[2*k*stride].imag,
             Wkn1.real, Wkn1.imag,
             x[2*k*stride + stride].real, x[2*k*stride + stride].imag);
      
      Xk1 = c_add(x[2*k*stride], c_mult(Wkn1, x[2*k*stride + stride]));

      complex Wkn2 = {cos(2*pi/N * (k + N/2)), -sin(2*pi/N * (k + N/2))};      
      printf("(%f,%f) + (%f,%f)*(%f,%f)\n",
             x[2*k*stride].real, x[2*k*stride].imag,
             Xk2.real, Xk2.imag,
             x[2*k*stride + stride].real, x[2*k*stride + stride].imag);
             
      Xk2 = c_add(x[2*k*stride], c_mult(Wkn2, x[2*k*stride + stride]));

      x[k*stride] = Xk1;
      printf("k*stride = %d\n", k*stride);
      x[(k + N/2)*stride] = Xk2;
      printf("(k + N/2)*stride = %d\n", (k + N/2)*stride);

    printf("[AFTER] Even FFT: ");
    for(int i=0; i<N/2; i++)
      printf("%d(%f,%f), ", i*stride*2, x[i*stride*2].real, x[i*stride*2].imag);
    printf("\n");

    printf("[AFTER] Odd FFT: ");
    for(int i=0; i<N/2; i++)
      printf("%d(%f,%f), ", i*stride*2 + stride, x[i*stride*2 + stride].real, x[i*stride*2 + stride].imag);
    printf("\n");

    }
  }
}

void fft_1d_slow(complex* x, int N, int stride) {
  complex y[N];
  for(int k=0; k<N; k++) {    
    y[k].real = 0;
    y[k].imag = 0;
    for(int n=0; n<N; n++) {
      complex Wkn = {cos(2*pi/N * k*n), -sin(2*pi/N * k*n)};
      y[k] = c_add(y[k], c_mult(x[n], Wkn));
    }
  }
  
  for(int k=0; k<N; k++)
    x[k] = y[k];
}
