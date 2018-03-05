#include <math.h>
#include "zetahyb.h"
#include <stdio.h>
#include "mpi.h"
#include "omp.h"
#include <stdlib.h>
#include <assert.h>


double zetahyb_function(int n, int mpi_size, int mpi_rank, int n_threads)
{
  // Number of iterations
  int iterations = n / mpi_size;

  // Initializing vectors
  double *vectors;

  // Calculating
  if(mpi_rank == 0){
      // Allocating space
      vectors = calloc((iterations * mpi_size) + mpi_size, sizeof(double));
      #pragma omp parallel for num_threads(n_threads)
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

int main(int argc, char *argv[]){

  int n = atoi(argv[2]);
  int n_threads = atoi(argv[1]);
  //init mpi
  int mpi_size, mpi_rank;
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  assert(ceil(log2(mpi_size)) == floor(log2(mpi_size)));

  double res = zetahyb_function(n, mpi_size, mpi_rank, n_threads);

  if(mpi_rank == 0){
    printf("Resiult: %f\n", res);
  }

  //Finalize the MPI
  MPI_Finalize();
}
