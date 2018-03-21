#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "mach1.h"
#include <math.h>
#include "mpi.h"
#define M_PI 3.14159265358979323846

int verification_mach1(int dist)
{
    int mpi_size, mpi_rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    assert(ceil(log2(mpi_size)) == floor(log2(mpi_size)));

    FILE *f = fopen("verification_results.txt", "w");
    double time1 = 0.0;
    if(mpi_rank == 0){
        time1 = MPI_Wtime();
    }
    for (int i = 1; i <= 24; i++){
        double mach = mach1_function(pow(2, i), mpi_size, mpi_rank, dist);

        if(mpi_rank == 0){
            printf("%f\n", mach);
            double error = fabs(M_PI - mach);
            fprintf(f, "N = %d - Error: %e\n", i, error);
        }
    }
    if(mpi_rank == 0){
        double time = MPI_Wtime() - time1;
        printf("Time: %e\n", time);
    }
    MPI_Finalize();

    return 0;
    }

int main(int argc, char *argv[])
{
    int dist = atoi(argv[1]);
    int ret = 0;
    ret |= verification_mach1(dist);
    return ret;
}