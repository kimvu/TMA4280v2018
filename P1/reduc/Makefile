
CC=mpicc
CFLAGS+=-std=c99 -O3
LDLIBS=-lm
TARGETS=zetareduc
NP=8

all: ${TARGETS}

clean:
	-rm -f ${TARGETS}

run:
	mpirun -np ${NP} zetareduc 4 0
