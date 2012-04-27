#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <string.h>
#include <math.h>
#include "virtual_network.h"
#include "sar.h"

//================================================================================
//======================== Nework Functions ======================================
//================================================================================

// Dummy row and column send operations
void send_message(int src_x, int src_y, int dest_x, int dest_y,
                  char* message, int size, char* send_buf){
  //Encode
  int header_size = 2*sizeof(int);
  int* tag = (int*)send_buf;
  tag[0] = src_x;
  tag[1] = src_y;
  //Copy Message
  memcpy(send_buf + header_size, message, size);
  //Send
  send_virtual_message(dest_x * Ny + dest_y, send_buf, size + header_size);
}

char* receive_message(int threadid, int* src_x, int* src_y, int* size){
  char* msg = receive_virtual_message(threadid);
  //Decode tag
  int header_size = 2*sizeof(int);
  int* tag = (int*)msg;
  *src_x = tag[0];
  *src_y = tag[1];
  *size = get_message_size(threadid) - header_size;
  return msg + header_size;
}

void free_message(int threadid){
  free_virtual_message(threadid);
}

//================================================================================
//========================= Send/Receive Row/Col =================================
//================================================================================

void send_row(int src_x, int src_y, char* message, int size, char* send_buf){
  if(src_x > 0)
    send_message(src_x, src_y, src_x - 1, src_y, message, size, send_buf);
  if(src_x < Nx-1)
    send_message(src_x, src_y, src_x + 1, src_y, message, size, send_buf);
}

void send_col(int src_x, int src_y, char* message, int size, char* send_buf){
  if(src_y > 0)
    send_message(src_x, src_y, src_x, src_y-1, message, size, send_buf);
  if(src_y < Ny-1)
    send_message(src_x, src_y, src_x, src_y+1, message, size, send_buf);
}

char* receive_and_forward(int ant_x, int ant_y, int* src_x, int* src_y, int* size,
                          int threadid, char* send_buf){
  char* msg = receive_message(threadid, src_x, src_y, size);
  //Forward Right
  if(*src_x < ant_x){
    if(ant_x < Nx - 1)
      send_message(*src_x, *src_y, ant_x + 1, ant_y, msg, *size, send_buf);
  }
  //Forward Left
  else if(*src_x > ant_x){
    if(ant_x > 0)
      send_message(*src_x, *src_y, ant_x - 1, ant_y, msg, *size, send_buf);
  }
  //Forward Up
  else if(*src_y < ant_y){
    if(ant_y < Ny - 1)
      send_message(*src_x, *src_y, ant_x, ant_y+1, msg, *size, send_buf);
  }
  //Forward Down
  else if(*src_y > ant_y){
    if(ant_y > 0)
      send_message(*src_x, *src_y, ant_x, ant_y-1, msg, *size, send_buf);
  }
  //Diagonal
  else {
    printf("[receive_and_forward] Diagonal message send not supported.\n");
    exit(-1);
  }
  return msg;
}

//================================================================================
//================================================================================
//================================================================================

void dump_data_buffer(char* data_buffer){
  int* id = (int*)data_buffer;
  complex* s = (complex*)(data_buffer + 2*sizeof(int));

  printf("id = %d, %d\n", id[0], id[1]);
  for(int i=0; i<Nf; i++){
    unsigned int* words = (unsigned int*)&s[i];
    printf("[%d, %d, %d]\n", words[0] >> 16, words[0] & 0xFFFF, words[1]);
  }
}

complex* Wkn_fft;
complex* Wkn_ifft;
complex* file_buffer;
complex* out_buffer;

void* sar_main(int threadid){
  // Get Position
  int ant_x = threadid / Ny;
  int ant_y = threadid % Ny;

  // Allocate space for send buffer
  char* send_buf = (char*)safe_malloc(MAX_MSG_SZ, "Failure to allocate send buffer.");

  // Allocate space 
  complex* s_local_1 = (complex*)safe_malloc(Nf * sizeof(complex), "Failed to initialize s.");
  complex* s_local_2 = (complex*)safe_malloc(Nf * sizeof(complex), "Failed to initialize s.");


  // Gather our data
  for(int i=0; i<Nf; i++){
    s_local_1[i] = file_buffer[(ant_x * 128/Nx) * 128 * 256 + (ant_y * 128/Ny) * 256 + (i * 256/Nf)];
  }
  
  // Generate some fake data
  /*
  for(int i=0; i<Nf; i++){
    int* words = (int*)(&s_local_1[i]);
    words[0] = (ant_x << 16) + ant_y;
    words[1] = i;
  }
  */

  // Send my local data to all processors in row
  // Every processor (x,y) sends an array, s_local_1, of size Nf, where
  // s_local_1[k] corresponds to s[x, y, k]
  send_row(ant_x, ant_y, (char*)s_local_1, Nf * sizeof(complex), send_buf);

  // Copy our own data to the correct location
  s_local_1[ant_x] = s_local_1[2*ant_x];
  s_local_1[ant_x + Nf/2] = s_local_1[2*ant_x + 1];
  
  // Receive neighbours data
  int step_1_count = 0;
  int step_2_count = 0;
  int step_3_count = 0;
  int step_4_count = 0;
  int step_5_count = 0;
  int step_6_count = 0;
  
  while(step_1_count < Nx - 1) {
    int src_x, src_y, size;
    complex* msg = (complex*)receive_and_forward(ant_x, ant_y, &src_x, &src_y, &size, threadid, send_buf);

    if(src_y == ant_y){
      //Step 1
      complex v1 = msg[2*ant_x];
      complex v2 = msg[2*ant_x + 1];
      s_local_1[src_x] = v1;
      s_local_1[src_x+Nf/2] = v2;
      step_1_count++;
    }
    else if(src_x == ant_x){
      //Step 2
      complex v1 = msg[ant_y];
      complex v2 = msg[ant_y+Nf/2];
      s_local_2[src_y] = v1;
      s_local_2[src_y + Nf/2] = v2;
      step_2_count++;
    }

    free_message(threadid);
  }

  // Do 2 1D FFTs (one for each frequency band) each of Nx points  
  fft_1d(s_local_1,      Nx, 1, Wkn_fft);
  fft_1d(s_local_1 + Nx, Nx, 1, Wkn_fft);

  // Print out some values
  /*
  if(ant_x == 7 && ant_y == 20){
    for(int i=0; i<Nf; i++){
      int* word = (int*)(&s_local_1[i]);
      printf("[x=%d, y=%d, f=%d]\n", word[0] >> 16, word[0] & 0xFFFF, word[1]);
    }
    printf("[Step 1] Data exchange across row complete\n");
  }
  */

  // Send my local data to all processors in column
  send_col(ant_x, ant_y, (char*)s_local_1, Nf * sizeof(complex), send_buf);

  // Copy out own data
  s_local_2[ant_y] = s_local_1[ant_y];
  s_local_2[ant_y + Nf/2] = s_local_1[ant_y + Nf/2];

  // Receive neighbours data
  while(step_2_count < Ny - 1) {
    int src_x, src_y, size;
    complex* msg = (complex*)receive_and_forward(ant_x, ant_y, &src_x, &src_y, &size, threadid, send_buf);
    if(src_x == ant_x){
      //Step 2
      complex v1 = msg[ant_y];
      complex v2 = msg[ant_y+Nf/2];
      s_local_2[src_y] = v1;
      s_local_2[src_y + Nf/2] = v2;
      step_2_count++;
    }
    else if(src_y == ant_y){
      //Step 3
      complex v1 = msg[ant_x];
      complex v2 = msg[ant_x + Nf/2];
      s_local_1[2*src_x] = v1;
      s_local_1[2*src_x+1] = v2;
      step_3_count++;
    }
    free_message(threadid);
  }

  // Print out some values
  /*
  if(ant_x == 7 && ant_y == 20){
    for(int i=0; i<Nf; i++){
      int* word = (int*)(&s_local_2[i]);
      printf("[x=%d, y=%d, f=%d]\n", word[0] >> 16, word[0] & 0xFFFF, word[1]);
    }
    printf("[Step 2] Data exchange across columns complete\n");
  }
  */

  // Do 2 1D FFT (one for each frequency band) each of Ny points
  fft_1d(s_local_2, Ny, 1, Wkn_fft);
  fft_1d(s_local_2+Ny, Ny, 1, Wkn_fft);
  
  // Send my local data to all processors in row
  send_row(ant_x, ant_y, (char*)s_local_2, Nf * sizeof(complex), send_buf);

  // Copy out own data
  s_local_1[2*ant_x] = s_local_2[ant_x];
  s_local_1[2*ant_x+1] = s_local_2[ant_x + Nf/2];

  // Receive neighbours data
  while(step_3_count < Nx - 1) {
    int src_x, src_y, size;
    complex* msg = (complex*)receive_and_forward(ant_x, ant_y, &src_x, &src_y, &size, threadid, send_buf);
    if(src_y == ant_y){
      //Step 3
      complex v1 = msg[ant_x];
      complex v2 = msg[ant_x + Nf/2];
      s_local_1[2*src_x] = v1;
      s_local_1[2*src_x+1] = v2;
      step_3_count++;
    }
    else if(src_x == ant_x){
      //Step 4
      complex v1 = msg[2*ant_y];
      complex v2 = msg[2*ant_y+1];
      s_local_2[src_y] = v1;
      s_local_2[src_y + Nf/2] = v2;
      step_4_count++;
    }
    free_message(threadid);
  }

  // Print out some values
  /*
  if(ant_x == 7 && ant_y == 20){
    for(int i=0; i<Nf; i++){
      int* word = (int*)(&s_local_1[i]);
      printf("[x=%d, y=%d, f=%d]\n", word[0] >> 16, word[0] & 0xFFFF, word[1]);
    }
    printf("[Step 3] Data exchange across rows complete\n");
  }
  */

  // Multiply each element in the frequency-domain signal by the
  // downward continuation phase operator.
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
  float kx = ant_y < Nx/2 ?
    2*pi/Dx * ant_y/Nx :
    2*pi/Dx * (ant_y - Nx)/Nx;
  float ky = ant_x < Ny/2 ?
    2*pi/Dy * ant_x/Ny :
    2*pi/Dy * (ant_x - Ny)/Ny;
  for(int n=0; n<Nf; n++){
    float w = 2*pi*(f0 + n*Df);
    float k = w/c_speed;
    float kz = sqrt(4*k*k - kx*kx - ky*ky);
    
    complex phi = c_jexp(kz * z0);
    s_local_1[n] = c_mult(s_local_1[n], phi);
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
  
  float n_interp[Nf];
  for(int n=0; n<Nf; n++){
    float kz = kz_min + (kz_max - kz_min) * n/(Nf - 1);
    float k = 0.5 * sqrt(kx*kx + ky*ky + kz*kz);
    n_interp[n] = (c_speed*k/(2*pi) - f0)/Df;
  }
  resample_1d(s_local_1, Nf, 1, n_interp);
  
  
  // Do a 1D IFFT over entire array of Nf points
  ifft_1d(s_local_1, Nf, 1, Wkn_ifft);
  
  send_col(ant_x, ant_y, (char*)s_local_1, Nf * sizeof(complex), send_buf);

  // Copy out own data
  s_local_2[ant_y] = s_local_1[2*ant_y];
  s_local_2[ant_y + Nf/2] = s_local_1[2*ant_y+1];

  // Receive neighbours data
  while(step_4_count < Ny - 1) {
    int src_x, src_y, size;
    complex* msg = (complex*)receive_and_forward(ant_x, ant_y, &src_x, &src_y, &size, threadid, send_buf);
    if(src_x == ant_x){
      //Step 4
      complex v1 = msg[2*ant_y];
      complex v2 = msg[2*ant_y+1];
      s_local_2[src_y] = v1;
      s_local_2[src_y + Nf/2] = v2;
      step_4_count++;
    }
    else if(src_y == ant_y){
      //Step 5
      complex v1 = msg[ant_x];
      complex v2 = msg[ant_x + Nf/2];
      s_local_1[src_x] = v1;
      s_local_1[src_x+Nf/2] = v2;
      step_5_count++;
    }
    free_message(threadid);
  }

  // Print out some values
  /*
  if(ant_x == 7 && ant_y == 20){
    for(int i=0; i<Nf; i++){
      int* word = (int*)(&s_local_2[i]);
      printf("[x=%d, y=%d, f=%d]\n", word[0] >> 16, word[0] & 0xFFFF, word[1]);
    }
    printf("[Step 4] Data exchange across columns complete\n");
  }
  */


  // Do 2 1D FFTs, each of Nx points
  ifft_1d(s_local_2, Nx, 1, Wkn_ifft);
  ifft_1d(s_local_2+Nx, Nx, 1, Wkn_ifft);

  send_row(ant_x, ant_y, (char*)s_local_2, Nf * sizeof(complex), send_buf);

  // Copy out own data
  s_local_1[ant_x] = s_local_2[ant_x];
  s_local_1[ant_x + Nf/2] = s_local_2[ant_x + Nf/2];

  // Receive neighbours data
  while(step_5_count < Nx - 1) {
    int src_x, src_y, size;
    complex* msg = (complex*)receive_and_forward(ant_x, ant_y, &src_x, &src_y, &size, threadid, send_buf);
    if(src_y == ant_y){
      //Step 5
      complex v1 = msg[ant_x];
      complex v2 = msg[ant_x + Nf/2];
      s_local_1[src_x] = v1;
      s_local_1[src_x+Nf/2] = v2;
      step_5_count++;
    }
    else if(src_x == ant_x){
      //Step 6
      complex v1 = msg[ant_y];
      complex v2 = msg[ant_y+Nf/2];
      s_local_2[2*src_y] = v1;
      s_local_2[2*src_y + 1] = v2;
      step_6_count++;
    }
    free_message(threadid);
  }  

  // Print out some values
  /*
  if(ant_x == 7 && ant_y == 20){
    for(int i=0; i<Nf; i++){
      int* word = (int*)(&s_local_1[i]);
      printf("[x=%d, y=%d, f=%d]\n", word[0] >> 16, word[0] & 0xFFFF, word[1]);
    }
    printf("[Step 5] Data exchange across rows complete\n");
  }
  */

  // Do 2 1D FFts, each of Ny points
  ifft_1d(s_local_1, Ny, 1, Wkn_ifft);
  ifft_1d(s_local_1+Ny, Ny, 1, Wkn_ifft);

  send_col(ant_x, ant_y, (char*)s_local_1, Nf * sizeof(complex), send_buf);

  // Copy out own data
  s_local_2[2*ant_y] = s_local_1[ant_y];
  s_local_2[2*ant_y + 1] = s_local_1[ant_y + Nf/2];

  // Receive neighbours data
  while(step_6_count < Ny - 1) {
    int src_x, src_y, size;
    complex* msg = (complex*)receive_and_forward(ant_x, ant_y, &src_x, &src_y, &size, threadid, send_buf);
    if(src_x == ant_x){
      //Step 6
      complex v1 = msg[ant_y];
      complex v2 = msg[ant_y+Nf/2];
      s_local_2[2*src_y] = v1;
      s_local_2[2*src_y + 1] = v2;
      step_6_count++;
    }
    else if(src_y == ant_y){
      //Step 7
      printf("Unreachable Statement\n");
    }
    free_message(threadid);
  }

  // Print out some values
  
  /* if(ant_x == 7 && ant_y == 20){ */
  /*   for(int i=0; i<Nf; i++){ */
  /*     c_print(s_local_2[i]); */
  /*     printf("\n"); */
  /*     /\* */
  /*     int* word = (int*)(&s_local_2[i]); */
  /*     printf("[x=%d, y=%d, f=%d]\n", word[0] >> 16, word[0] & 0xFFFF, word[1]); */
  /*     *\/ */
  /*   } */
  /*   printf("[Step 6] Data exchange across columns complete\n"); */
  /* } */
  

  // Write out our data
  for(int i=0; i<Nf; i++)
    out_buffer[ant_x * Ny * Nf + ant_y * Nf + i] = s_local_2[i];
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
        data[n] = num;
        n++;
      }

  fclose(file);
}

int main(){
  // Read in file
  file_buffer = (complex*)safe_malloc(256 * 128 * 128 * sizeof(complex), "Failed to allocate file buffer.\n");
  out_buffer = (complex*)safe_malloc(Nx * Ny * Nf * sizeof(complex), "Failed to allocate out buffer.\n");
  read_original_data(file_buffer, "scene_4.dat");

  // Compute FFT coefficients
  Wkn_fft = precompute_fft_coefficients();
  Wkn_ifft = precompute_ifft_coefficients();
  
  init_virtual_network(Nx, Ny);
  start_virtual_network(&sar_main);
  free_virtual_network();

  // Write out file
  write_data(out_buffer, "network_sar.out");
  return 0;
}
