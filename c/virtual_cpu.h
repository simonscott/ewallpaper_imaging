#ifndef VIRTUAL_CPU_H
#define VIRTUAL_CPU_H
#include <pthread.h>

//================================================================================
//============================ CPU Ports =========================================
//================================================================================

typedef struct {
  // Message Buffer & Status
  char* buffer;
  int space_available;
  int message_waiting;
  // Notifier
  pthread_cond_t signal;
  pthread_mutex_t signal_lock;
} cpu_port;

void init_cpu_port(cpu_port* p);
void notify_port_listeners(cpu_port* port);
void send_port(cpu_port* port, void* message, int message_size);
char* receive_port(cpu_port* port);
void free_port(cpu_port* port, void* buffer);

//================================================================================
//=============================== CPU ============================================
//================================================================================

typedef struct {
  // Message Ports
  cpu_port left_port;
  cpu_port right_port;
  cpu_port up_port;
  cpu_port down_port;
  // Execution Thread
  pthread_t thread;
} cpu;

// Directions
extern const int up;
extern const int down;
extern const int left;
extern const int right;

void init_cpu(cpu* processor);
void start_cpu(cpu* processor, void* (*cpu_thread)(void*), void* args);
cpu_port* get_port(cpu* processor, int direction);

//================================================================================
//============================ Directions ========================================
//================================================================================
int opposite_dir(int dir);
void step_in_dir(int dir, int* x, int* y);

#endif
