
# This Makefile is still simple and the compilations lines are just duplicated.
# This is suboptimal but good enough for this case.
CC=gcc
TARGET=mach0
SRC=mach0.c

# Notice that the prerequisites can be stored in a variable.
all:	$(SRC)
	$(CC) -c mach0.c
	$(CC) -c mach0test.c
	$(CC) -o $(TARGET)test mach0test.o mach0.o -lm

clean:
	rm -f *.o $(TARGET)

utest: mach0test
	./mach0test
