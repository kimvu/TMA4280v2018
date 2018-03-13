CC=mpicc
TARGET=mach1
SRC=mach1.c

all:	$(SRC)
	$(CC) -c mach1.c
	$(CC) -c mach1verification.c
	$(CC) -o $(TARGET)verification mach1verification.o  mach1.o -lm

clean:
	rm -f *.o $(TARGET)

vtest: mach1verification
	mpirun -np 16 ./mach1verification 0
