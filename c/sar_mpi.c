#include "sar.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//================================================================================
//==================== Network Functions =========================================
//================================================================================

// Sends the provided message to all processors in the same row as this processor.
// This is implemented using MPI Broadcast.
void send_row(void* message, int size)
{
}

// Sends the provided message to all processors in the same column as this processor.
// This is implemented using MPI Broadcast.
void send_col(void* message, int size)
{
}

// Receives a message from any processor in the same row as this processor.
// The message is placed in the provided buffer, and the size of the message (in bytes)
// is returned.
// This is implemented using an MPI Blocking Read.
int receive_row(void* buffer)
{
  return 0;
}

// Receives a message from any processor in the same column as this processor.
// The message is placed in the provided buffer, and the size of the message (in bytes)
// is returned.
// This is implemented using an MPI Blocking Read.
int receive_col(void* buffer)
{
  return 0;
}

//================================================================================
//==================== The Main Function =========================================
//================================================================================

int main(int argc, char* argv[])
{
  // Declare local variables

  // MPI setup

  // Determine position in 2D processor array

  // 
}
