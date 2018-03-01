#include <math.h>
#include "zeta1.h"
#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>


double zeta1_function(int n, int mpi_size, int mpi_rank)
{
  // Number of iterations
  int iterations = n / mpi_size + 1;

  // Initializing vectors
  double *vectors;

  // Calculating
  if(mpi_rank == 0){
      // Allocating space
      vectors = calloc((iterations * mpi_size) + mpi_size, sizeof(double));
      for (int i=1; i<=n; i++) {
          vectors[i-1] = 1.0/((double)i*(double)i);
      }
  }
  double *local_values = calloc(iterations, sizeof(double));
  // Partitions the array
  MPI_Scatter(vectors, iterations, MPI_DOUBLE, local_values, iterations, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Calculating local sum, reducing to mpi_rank 0
  double values_sum = 0.0;
  for (int i=0; i<iterations; i++) {
      values_sum += local_values[i];
  }

  double sum = 0.0;
  MPI_Reduce(&values_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  // Freeing memory
  if(mpi_rank == 0){
      free (vectors);
  }
  free(local_values);

  return sqrt(sum*6);
}
