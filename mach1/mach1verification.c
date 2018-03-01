#include <stdio.h>
#include <assert.h>
#include "mach1.h"
#include <math.h>
#include "mpi.h"
#define M_PI 3.14159265358979323846

int verification_mach1()
{
    int mpi_size, mpi_rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    FILE *f = fopen("results/verification_results.txt", "w");
    double errors[24];
    for (int i = 1; i <= 24; i++){
        double time1 = MPI_Wtime();
        double mach = mach1_function(pow(2, i), mpi_size, mpi_rank);
        double time2 = MPI_Wtime();

        double time = time2 - time1;
        double error = (fabs(M_PI - mach));
        fprintf(f, "Error: %e - Time: %e\n", error, time);
        printf("Error: %e - Time: %e\n", error, time);
    }

    MPI_Finalize();

    return 0;
    }

int main(int argc, char *argv[])
{
    int ret = 0;
    ret |= verification_mach1();
    return ret;
}
