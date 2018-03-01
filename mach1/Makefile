
# This Makefile is still simple and the compilations lines are just duplicated.
# This is suboptimal but good enough for this case.
CC=mpicc
TARGET=mach1
SRC=mach1.c

# Notice that the prerequisites can be stored in a variable.
all:	$(SRC)
	$(CC) -c mach1.c
	$(CC) -c mach1verification.c
	$(CC) -o $(TARGET)verification mach1verification.o  mach1.o -lm

clean:
	rm -f *.o $(TARGET)

vtest: mach1verification
	mpirun -np 4 ./mach1verification
