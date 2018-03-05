
# This Makefile is still simple and the compilations lines are just duplicated.
# This is suboptimal but good enough for this case.
CC=mpicc
TARGET=zetahyb
SRC=zetahyb.c
LDLIBS=-lm
# Notice that the prerequisites can be stored in a variable.
all:	$(TARGET)

clean:
	rm -f *.o $(TARGET)

run:
	mpirun -np 4 zetahyb 16 1024
