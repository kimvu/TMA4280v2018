#include <math.h>
#include "zetareduc.h"
#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#include <assert.h>


double zetareduc_function(int n, int mpi_size, int mpi_rank, int dist)
{

  // Number of iterations
  int iterations = n / mpi_size;

  // Initializing vectors
  double *vectors;
  double the_pi = 0.0;
  double values_sum = 0.0;
  // Calculating
  if(dist)
  {
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
  the_pi = recursive_allreduce(mpi_size, mpi_rank, values_sum);
  return sqrt(the_pi*6);
}


// Global reduction
double allreduce(int mpi_size, int mpi_rank, double values_sum){
    double sum = values_sum;
    MPI_Allreduce(&values_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sum;
}

// Recursive
double recursive_allreduce(int mpi_size, int mpi_rank, double values_sum){
    double sum = values_sum;
    for(int i = 1; i < mpi_size; i = i << 1){
        int exchange_value = mpi_rank ^ i;

        double recv = 0;
        MPI_Sendrecv(&values_sum, 1, MPI_DOUBLE, exchange_value, 0,
                     &recv, 1, MPI_DOUBLE, exchange_value, MPI_ANY_TAG,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        sum += recv;
    }
    return sum;

}

int main(int argc, char **argv){

    // Initializing MPI
    int mpi_size, mpi_rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    // Getting argument
    int n = atoi(argv[1]);
    int dist = atoi(argv[2]);

    double sum = zetareduc_function(n, mpi_size, mpi_rank, dist);
    printf("MPI Rank: %d - Result: %f\n", mpi_rank, sum);

    // Finalizing MPI
    MPI_Finalize();

    exit(EXIT_SUCCESS);

}
