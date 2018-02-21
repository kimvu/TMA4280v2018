
# This Makefile is still simple and the compilations lines are just duplicated.
# This is suboptimal but good enough for this case.
CC=gcc
TARGET=zeta0
SRC=zeta0.c

# Notice that the prerequisites can be stored in a variable.
all:	$(SRC)
	$(CC) -c zeta0.c
	$(CC) -c zeta0test.c
	$(CC) -c zeta0verification.c
	$(CC) -o $(TARGET)test zeta0test.o zeta0.o -lm
	$(CC) -o $(TARGET)verification zeta0verification.o  zeta0.o -lm

clean:
	rm -f *.o $(TARGET)

utest: zeta0test
	./zeta0test

vtest: zeta0verification
	./zeta0verification
