#include "mach1.h"
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>

void arctan(int n, double *vectors)
{
  double sum = 0;

  for (int i = 1; i <= n; i++){
    double term1 = pow(-1, i-1)*(pow((double)0.2, 2*i-1)/(2*i-1));
    double term2 = pow(-1, i-1)*(pow((double)1/239, 2*i-1)/(2*i-1));
    vectors[i] = 16*term1 - 4*term2;
  }
}

double mach1_function(int n, int mpi_size, int mpi_rank, int dist){

  // Number of iterations
  int iterations = n / mpi_size;
  double values_sum = 0.0;
  // Initializing vectors
  double *vectors;
  if (dist)
  {
    // Calculating
    if(mpi_rank == 0){
        // Allocating space
        vectors = calloc((iterations * mpi_size) + mpi_size, sizeof(double));
        arctan(n,vectors);
    }

    double *local_values = calloc(iterations, sizeof(double));
    // Partitions the array
    MPI_Scatter(vectors, iterations, MPI_DOUBLE, local_values, iterations, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Calculating local sum, reducing to mpi_rank 0

    for (int i=0; i<iterations; i++) {
        values_sum += local_values[i];

    }

    // Freeing memory
    if(mpi_rank == 0){
        free (vectors);
    }
    free(local_values);
  }
  else
  {
    // Calculating local sum, reducing to mpi_rank 0
    int i0 = iterations * mpi_rank;
    int i1 = i0 + iterations;
    for (int i=i0; i<i1; i++) {
      double term1 = pow(-1, i-1)*(pow((double)0.2, 2*i-1)/(2*i-1));
      double term2 = pow(-1, i-1)*(pow((double)1/239, 2*i-1)/(2*i-1));
      values_sum += 16*term1 - 4*term2;
    }
  }

  double the_pi = 0.0;
  MPI_Reduce(&values_sum, &the_pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  return the_pi;
}
