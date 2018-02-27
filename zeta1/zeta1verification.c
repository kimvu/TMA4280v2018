#include <stdio.h>
#include <assert.h>
#include "zeta1.h"
#include <math.h>
#include <mpi.h>
#define M_PI 3.14159265358979323846

int verification_zeta1()
{
  int mpi_size, mpi_rank;
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  double errors[24];
  for (int i = 1; i <= 24; i++){
    errors[i-1] = (fabs(M_PI - zeta1_function(pow(2, i), mpi_size, mpi_rank)));
  }

  FILE *f = fopen("results/verification_results.txt", "w");

  for (int i = 1; i <= 24; i++){
    printf("%e\n", errors[i-1]);
    fprintf(f, "%e\n", errors[i-1]);
  }
  MPI_Finalize();

  return 0;
}

int main(int argc, char *argv[])
{
    int ret = 0;
    ret |= verification_zeta1();
    return ret;
}
