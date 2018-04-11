#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

// Send counts and displacements
int *sendcounts;
int *recvcounts;
int *senddispl;
int *recvdispl;

MPI_Datatype matrixcolumn;
MPI_Datatype matrixcolumntype;

MPI_Datatype sendtypelarge;
MPI_Datatype sendtypesmall;
MPI_Datatype recvtypelarge;
MPI_Datatype recvtypesmall;