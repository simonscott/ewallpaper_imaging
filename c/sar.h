#ifndef SAR_H
#define SAR_H

//Network Functions
void send_row(void* message, int size);
void send_col(void* message, int size);
void* receive_row();
void* receive_col();

//Complex Types
typedef struct {
  float real;
  float imag;
} complex;

//Complex Number Functions
complex c_mult(complex x, complex y);
complex c_add(complex x, complex y);

//Signal Processing Functions
void fft_1d(complex* x, int n, int stride);
void fft_1d_slow(complex* x, int n, int stride);
void ifft_1d(complex* x, int n, int stride);

//Application Specific Math
//void interp_1d(complex* x, int n, ?);
void precompute_phase_shifts(complex* x, int n);
void element_wise_mult(complex* x, complex* y, int n);

//Simulation setup and teardown
void init_simulation(complex* data, char* filename);
void end_simulation();

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

#endif
