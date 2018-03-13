
CC=gcc
TARGET=mach0
SRC=mach0.c

all:	$(SRC)
	$(CC) -c mach0.c
	$(CC) -c mach0test.c
	$(CC) -c mach0verification.c
	$(CC) -o $(TARGET)test mach0test.o mach0.o -lm
	$(CC) -o $(TARGET)verification mach0verification.o  mach0.o -lm

clean:
	rm -f *.o $(TARGET)

utest: mach0test
	./mach0test

vtest: mach0verification
	./mach0verification
