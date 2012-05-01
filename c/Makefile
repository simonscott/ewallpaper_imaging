CC = CC
MPCC = CC
CFLAGS = -O3 -c99
#CFLAGS = -g
LIBS = -lm

TARGETS = sar_serial sar_mpi new_mpi

all: $(TARGETS)

sar_serial: sar.o sar_serial.o
	$(CC) -o $@ sar.o sar_serial.o $(LIBS)
sar_mpi: sar.o sar_mpi.o virtual_network.o
	$(MPCC) -o $@ -lpthread $(LIBS) $(MPILIBS) -o $@ sar.o sar_mpi.o virtual_network.o
new_mpi: virtual_cpu.o mpi_network.o simple_simulation.o
	$(MPCC) -o $@ -lpthread $(LIBS) $(MPILIBS) virtual_cpu.o mpi_network.o simple_simulation.o

sar_serial.o: sar_serial.c sar.c sar.h
	$(CC) -c $(CFLAGS) sar_serial.c
sar_mpi.o: sar_mpi.c sar.c sar.h virtual_network.c virtual_network.h
	$(MPCC) -c $(CFLAGS) sar_mpi.c
sar.o: sar.c sar.h
	$(CC) -c $(CFLAGS) sar.c
virtual_network.o: virtual_network.c virtual_network.h
	$(CC) -c $(CFLAGS) -lpthread virtual_network.c

virtual_cpu.o: virtual_cpu.c virtual_cpu.h
	$(CC) -c $(CFLAGS) -lpthread virtual_cpu.c
mpi_network.o: mpi_network.c virtual_cpu.h
	$(MPCC) -c $(CFLAGS) mpi_network.c
simple_simulation.o: simple_simulation.c virtual_cpu.h sar.h
	$(MPCC) -c $(CFLAGS) simple_simulation.c

clean:
	rm -f *.o $(TARGETS)




