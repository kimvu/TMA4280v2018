#include <math.h>
#include "zetahyb.h"
#include <stdio.h>
#include "mpi.h"
#include "omp.h"
#include <stdlib.h>
#include <assert.h>

//Hybrid function
double zetahyb_function(int n, int mpi_size, int mpi_rank, int n_threads, int dist)
{

  // Number of iterations
  int iterations = n / mpi_size;

  // Initializing vectors ++
  double *vectors;
  double the_pi = 0.0;
  double values_sum = 0.0;

  // Check if the work load should be distributed
  if(!dist)
  {
    if(mpi_rank == 0){
        // Allocating space
        vectors = calloc((iterations * mpi_size) + mpi_size, sizeof(double));

        // OpemMP
        #pragma omp parallel for num_threads(n_threads)
        for (int i=1; i<=n; i++) {
            vectors[i-1] = 1.0/((double)i*(double)i);
        }
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
    #pragma omp parallel for num_threads(n_threads)
    for (int i=i0; i<i1; i++) {
      values_sum += 1.0/((double)i*(double)i);
    }
  }

  // Collects the calculations from all the processes, and reduces to the variable the_pi
  MPI_Reduce(&values_sum, &the_pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  return sqrt(the_pi*6);
}

// The main function
int main(int argc, char *argv[]){

  int n = atoi(argv[2]);
  int n_threads = atoi(argv[1]);
  int dist = atoi(argv[3]);

  //initializing mpi
  int mpi_size, mpi_rank;
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  // Checks that the number of processes is a power of 2
  assert(ceil(log2(mpi_size)) == floor(log2(mpi_size)));

  double res = zetahyb_function(n, mpi_size, mpi_rank, n_threads, dist);

  if(mpi_rank == 0){
    printf("Result: %f\n", res);
  }

  //Finalize the MPI
  MPI_Finalize();
}
