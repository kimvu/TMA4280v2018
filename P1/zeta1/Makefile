
# This Makefile is still simple and the compilations lines are just duplicated.
# This is suboptimal but good enough for this case.
CC=mpicc
TARGET=zeta1
SRC=zeta1.c

# Notice that the prerequisites can be stored in a variable.
all:	$(SRC)
	$(CC) -c zeta1.c
	$(CC) -c zeta1verification.c
	$(CC) -o $(TARGET)verification zeta1verification.o  zeta1.o -lm

clean:
	rm -f *.o $(TARGET)

vtest: zeta1verification
	mpirun -np 16 ./zeta1verification
