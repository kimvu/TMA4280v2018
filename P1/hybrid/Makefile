CC=mpicc
TARGET=zetahyb
SRC=zetahyb.c
LDLIBS=-lm

all:	$(TARGET)

clean:
	rm -f *.o $(TARGET)

run:
	mpirun -np 8 zetahyb 100 10 0
