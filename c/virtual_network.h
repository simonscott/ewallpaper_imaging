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
// sends are nonblocking, receives are blocking.

// send_message :
// threadid is the id of the destination processor
// message is a pointer to the message to be sent
// size is the length of the message in bytes
// the message is copied from the message buffer to the network buffer.
void send_virtual_message(int threadid, char* message, int size);

// Receiving Messages
// receive_message must always be followed by free_message after the message has been
// processed. Otherwise the network will overflow.

// receive_message :
// Should always be called like this receive_message(MYTHREAD)
// Returns a pointer to the first message in the network buffer.
// If there is no message in the network buffer, then block until there is one.
char* receive_virtual_message(int threadid);

// free_message :
// Should always be called like this free_message(MYTHREAD, pointer_to_first_message)
// Frees the space allocated for the first message in the network buffer.
void free_virtual_message(int threadid, char* message);

#endif
