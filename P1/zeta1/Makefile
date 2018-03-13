
CC=mpicc
TARGET=zeta1
SRC=zeta1.c

all:	$(SRC)
	$(CC) -c zeta1.c
	$(CC) -c zeta1verification.c
	$(CC) -o $(TARGET)verification zeta1verification.o  zeta1.o -lm

clean:
	rm -f *.o $(TARGET)

vtest: zeta1verification
	mpirun -np 16 ./zeta1verification 0
