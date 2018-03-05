
CC=mpicc
CFLAGS+=-std=c99 -O3
LDLIBS=-lm
TARGETS=zetareduc
NP=16

all: ${TARGETS}

clean:
	-rm -f ${TARGETS}

run:
	mpirun -np ${NP} zetareduc 42