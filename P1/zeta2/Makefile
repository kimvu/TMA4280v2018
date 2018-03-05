
# This Makefile is still simple and the compilations lines are just duplicated.
# This is suboptimal but good enough for this case.
TARGET=zeta2
SRC=zeta2.c
LDLIBS=-lm
# Notice that the prerequisites can be stored in a variable.
all:	$(TARGET)

clean:
	rm -f *.o $(TARGET)

run:
	./zeta2 4 1024
