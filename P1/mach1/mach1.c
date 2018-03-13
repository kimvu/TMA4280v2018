#include "mach1.h"
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>

// The arctan function
void arctan(int n, double *vectors)
{
  double sum = 0;

  for (int i = 1; i <= n; i++){
    double term1 = pow(-1, i-1)*(pow((double)0.2, 2*i-1)/(2*i-1));
    double term2 = pow(-1, i-1)*(pow((double)1/239, 2*i-1)/(2*i-1));
    vectors[i] = 16*term1 - 4*term2;
  }
}

// The mach funtion
double mach1_function(int n, int mpi_size, int mpi_rank, int dist){

  // Number of iterations
  int iterations = n / mpi_size;
  double values_sum = 0.0;

  // Initializing vectors
  double *vectors;

  // Check if the work load should be distributed
  if (dist)
  {
    if(mpi_rank == 0){
        // Allocating space and calculating
        vectors = calloc((iterations * mpi_size) + mpi_size, sizeof(double));
        arctan(n,vectors);
    }

    double *local_values = calloc(iterations, sizeof(double));

    // Splitting data to different processes from the vector
    MPI_Scatter(vectors, iterations, MPI_DOUBLE, local_values, iterations, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Calculating sum
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

    // Arctan function for each process
    for (int i=i0; i<i1; i++) {
      double term1 = pow(-1, i-1)*(pow((double)0.2, 2*i-1)/(2*i-1));
      double term2 = pow(-1, i-1)*(pow((double)1/239, 2*i-1)/(2*i-1));
      values_sum += 16*term1 - 4*term2;
    }
  }

  double the_pi = 0.0;

  // Collects the calculations from all the processes, and reduces to the variable the_pi
  MPI_Reduce(&values_sum, &the_pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  return the_pi;
}
