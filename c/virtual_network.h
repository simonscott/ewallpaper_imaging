#ifndef VIRTUAL_NETWORK_H
#define VIRTUAL_NETWORK_H

// Starting a virtual network
// e.g. start_virtual_network(10, 10, &my_main)
// will start a virtual network with 10 by 10 processors on a grid, and initialize
// each processor to run my_main.
// my_main must take a single integer argument MYTHREAD, and return a void pointer.
typedef void* (*processor_main_function)(int MYTHREAD);
void start_virtual_network(int rows, int cols, processor_main_function p_main);

// Message Passing Functions
void send_message(int threadid, void* message, int size);
void* receive_message(int threadid);
void free_message(int threadid, void* message);

#endif
