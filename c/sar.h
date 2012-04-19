#ifndef SAR_H
#define SAR_H

//Physical Constants
#define c_speed 299792458
#define pi M_PI

//Antenna Physical Parameters
#define Nx 128
#define Ny 128
#define Dx 0.012
#define Dy 0.012

//Antenna RF Parameters
#define f0 10e9
#define B 2e9
#define Nf 256
#define Df (B/Nf)

//Scene Parameters
#define z0 1

//Complex Types
typedef struct {
  float real;
  float imag;
} complex;

//Network Functions
void send_row(void* message, int size);
void send_col(void* message, int size);
int receive_row(void* buffer);
int receive_col(void* buffer);

//Complex Number Functions
complex c_mult(complex x, complex y);
complex c_add(complex x, complex y);
complex c_jexp(float theta);
complex c_scalar_mult(complex x, float a);
complex c_scalar_div(complex x, float a);
void c_print(complex x);

//Signal Processing Functions
void fft_1d(complex* x, int N, int stride);
void ifft_1d(complex* x, int N, int stride);
void resample_1d(complex* x, int N, int stride, float* n);

//Application Specific Math
void precompute_phase_shifts(complex* x, int n);
void element_wise_mult(complex* x, complex* y, int n);

//Simulation setup and teardown
void read_data(complex* data, char* filename);
void write_data(complex* data, char* filename);

//Memory allocation and deallocation
void* safe_malloc(int size, char* error_message);

#endif
