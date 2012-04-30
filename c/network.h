#ifndef NETWORK_H
#define NETWORK_H

void send_message(int x, int y, int direction, void* message, int message_size);
char* receive_message(int x, int y, int direction);
void free_network_port(int x, int y, int direction, void* buffer);

#endif
