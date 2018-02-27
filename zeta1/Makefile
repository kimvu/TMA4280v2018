
# This Makefile is still simple and the compilations lines are just duplicated.
# This is suboptimal but good enough for this case.
CC=gcc
TARGET=zeta1
SRC=zeta1.c

# Notice that the prerequisites can be stored in a variable.
all:	$(SRC)
	$(CC) -c zeta1.c zeta1.o

clean:
	rm -f *.o $(TARGET)
