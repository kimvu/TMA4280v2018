#include "mach0.h"

#include <stdio.h>
#include <math.h>

double arctan(int n, double x)
{
  double sum = 0;

  for (double i = 1; i <= n; i++){
    sum += pow(-1, i-1)*(pow(x, 2*i-1)/(2*i-1));
  }

  return sum;
}

double mach_function(int n, int mpi_size, int mpi_rank){

  // Number of iterations
  int iterations = n / mpi_size + 1;

  // Initializing vectors
  double *vectors;

  // Calculating
  if(mpi_rank == 0){
      // Allocating space
      vectors = calloc((iterations * mpi_size) + mpi_size, sizeof(double));
      for (int i=1; i<=n; i++) {
        double first_term = 16*arctan(n, (double)0.2);
        double second_term = 4*arctan(n, (double)1/239);

        vectors[i-1] = first_term - second_term;
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

  double the_pi = 0.0;
  MPI_Reduce(&values_sum, &the_pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  // Freeing memory
  if(mpi_rank == 0){
      free (vectors);
  }
  free(local_values);

  return sum;
}
