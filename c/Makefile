CC = CC
MPCC = CC
CFLAGS = -O3 -c99
#CFLAGS = -g
LIBS = -lm

TARGETS = sar_serial sar_mpi

all: $(TARGETS)

sar_serial: sar.o sar_serial.o
	$(CC) -o $@ sar.o sar_serial.o $(LIBS)
sar_mpi: sar.o sar_mpi.o virtual_network.o
	$(MPCC) -o $@ -lpthread $(LIBS) $(MPILIBS) -o $@ sar.o sar_mpi.o virtual_network.o

sar_serial.o: sar_serial.c sar.c sar.h
	$(CC) -c $(CFLAGS) sar_serial.c
sar_mpi.o: sar_mpi.c sar.c sar.h virtual_network.c virtual_network.h
	$(MPCC) -c $(CFLAGS) sar_mpi.c
sar.o: sar.c sar.h
	$(CC) -c $(CFLAGS) sar.c
virtual_network.o: virtual_network.c virtual_network.h
	$(CC) -c $(CFLAGS) -lpthread virtual_network.c

clean:
	rm -f *.o $(TARGETS)

