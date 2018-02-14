
# This Makefile is still simple and the compilations lines are just duplicated.
# This is suboptimal but good enough for this case.
CC=gcc
TARGET=zeta0
SRC=zeta0.c

# Notice that the prerequisites can be stored in a variable.
all:	$(SRC)
	$(CC) -c zeta0.c
	$(CC) -o $(TARGET) zeta0.o

clean:
	rm -f *.o $(TARGET)
