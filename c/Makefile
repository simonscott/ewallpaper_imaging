CC = gcc
CFLAGS = -O3
LIBS = -lm

TARGETS = sar_serial sar_mpi

all: $(TARGETS)

sar_serial: sar.o sar_serial.o
	$(CC) -o $@ sar.o sar_serial.o $(LIBS)
sar_mpi: sar.o sar_mpi.o
	$(CC) -o $@ sar.o sar_mpi.o $(LIBS)

sar_serial.o: sar_serial.c sar.c sar.h
	$(CC) -c $(OPENMP) $(CFLAGS) sar_serial.c
sar_mpi.o: sar_mpi.c sar.c sar.h
	$(CC) -c $(OPENMP) $(CFLAGS) sar_mpi.c
sar.o: sar.c sar.h
	$(CC) -c $(CFLAGS) sar.c

clean:
	rm -f *.o $(TARGETS)

