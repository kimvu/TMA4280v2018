#include <math.h>
#include "zeta1.h"
#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>

// The zeta function
double zeta1_function(int n, int mpi_size, int mpi_rank, int dist)
{

  // Number of iterations
  int iterations = n / mpi_size;

  // Initializing vectors ++
  double *vectors;
  double the_pi = 0.0;
  double values_sum = 0.0;

    // Check if the work load should be distributed
  if(dist)
  {
    if(mpi_rank == 0){
        // Allocating space and calculating
        vectors = calloc((iterations * mpi_size) + mpi_size, sizeof(double));
        for (int i=1; i<=n; i++) {
            vectors[i-1] = 1.0/((double)i*(double)i);
        }
    }
    double *local_values = calloc(iterations, sizeof(double));

    // Splitting data to different processes from the vector
    MPI_Scatter(vectors, iterations, MPI_DOUBLE, local_values, iterations, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Calculating local sum
    for (int i=0; i<iterations; i++) {
        values_sum += local_values[i];
    }

    // Freeing memory
    if(mpi_rank == 0){
        free (vectors);
    }
    free(local_values);
  }else{
    int i0;
    int i1;
    if(mpi_size > n){ // if mpi_size are bigger than n
        if(mpi_rank < n){ // share the rest of the tasks
            i0 = mpi_rank + 1;
            i1 = i0 + 1;
        }else{ // do nothing
          i0 = 0;
          i1 = 0;
        }
    }else{
      i0 = iterations * mpi_rank + 1;
      i1 = iterations + i0 ;
    }
    for (int i=i0; i<i1; i++) {
      values_sum += 1.0/((double)i*(double)i);
    }
  }

  // Collects the calculations from all the processes, and reduces to the variable the_pi
  MPI_Reduce(&values_sum, &the_pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  return sqrt(the_pi*6);
}
