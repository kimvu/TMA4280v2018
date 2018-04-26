/**
 * C program to solve the two-dimensional Poisson equation on
 * a unit square using one-dimensional eigenvalue decompositions
 * and fast sine transforms.
 *
 * Einar M. Rønquist
 * NTNU, October 2000
 * Revised, October 2001
 * Revised by Eivind Fonn, February 2015
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define PI 3.14159265358979323846
#define true 1
#define false 0

typedef double real;
typedef int bool;

// Function prototypes
real *mk_1D_array(size_t n, bool zero);
real **mk_2D_array(size_t n1, size_t n2, bool zero);
void transpose(real **bt, real **b, size_t m);
real rhs(real x, real y);
real answer(real x, real y);

// Functions implemented in FORTRAN in fst.f and called from C.
// The trailing underscore comes from a convention for symbol names, called name
// mangling: if can differ with compilers.
void fst_(real *v, int *n, real *w, int *nn);
void fstinv_(real *v, int *n, real *w, int *nn);

// Send/receice counts and displacements
int *sendcounts;
int *recvcounts;
int *rdispls;
int *sdispls;

// from/to for each processor
int *from_proc;
int *to_proc;

// MPI Datatypes
MPI_Datatype m_column;
MPI_Datatype m_columntype;


int main(int argc, char **argv)
{
    // Initializing MPI
    int mpi_size, mpi_rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    if (argc < 2) {
        printf("Usage:\n");
        printf("  poisson n\n\n");
        printf("Arguments:\n");
        printf("  n: the problem size (must be a power of 2)\n");
    }

    // Start Time
    double time1;
    if(mpi_rank == 0){
        time1 = MPI_Wtime();
    }

    /*
     *  The equation is solved on a 2D structured grid and homogeneous Dirichlet
     *  conditions are applied on the boundary:
     *  - the number of grid points in each direction is n+1,
     *  - the number of degrees of freedom in each direction is m = n-1,
     *  - the mesh size is constant h = 1/n.
     */
    int n = atoi(argv[1]);
    int n_threads = atoi(argv[2]);
    int m = n - 1;
    real h = 1.0 / n;

    // Allocating memory, for load balancing
    from_proc = (int*) malloc (mpi_size * sizeof (int));
    to_proc = (int*) malloc (mpi_size * sizeof (int));

    // Allocating memory, for the transposing
    sendcounts = (int *)malloc( mpi_size * sizeof(int) );
    recvcounts = (int *)malloc( mpi_size * sizeof(int) );
    rdispls = (int *)malloc( mpi_size * sizeof(int) );
    sdispls = (int *)malloc( mpi_size * sizeof(int) );

    // Rows in each process
    int num_row = m / mpi_size;
    int row_rem = m % mpi_size;

    // Calculating from and to
    for (int i = 0; i < mpi_size; i++)
    {
      from_proc[i] = i * num_row + (i < row_rem ? i : row_rem);
      to_proc[i] = from_proc[i] + num_row + (i < row_rem ? 1 : 0);
    }

    //TODO HUSK Å SETTE DEM I FOR LØKKER
    int from = from_proc[mpi_rank];
    int to = to_proc[mpi_rank];

    // Assigns work load for each process
    for (int i = 0; i < mpi_size;i++){
      int num_row = to - from;
      sendcounts[i] = num_row * m;
      sdispls[i] = from * m;
      recvcounts[i] = to_proc[i] - from_proc[i];
      rdispls[i] = from_proc[i];

    // Initializing MPI matrix
      MPI_Type_vector (m, 1, m, MPI_DOUBLE, &m_column);
      MPI_Type_commit (&m_column);
      MPI_Type_create_resized (m_column, 0, sizeof(double), &m_columntype);
      MPI_Type_commit (&m_columntype);
    }

    /*
     * Grid points are generated with constant mesh size on both x- and y-axis.
     */
    real *grid = mk_1D_array(n+1, false);

#pragma omp parallel for num_threads(n_threads) schedule(static)
    for (size_t i = 0; i < n+1; i++) {
        grid[i] = i * h;
    }

    /*
     * The diagonal of the eigenvalue matrix of T is set with the eigenvalues
     * defined Chapter 9. page 93 of the Lecture Notes.
     * Note that the indexing starts from zero here, thus i+1.
     */
    real *diag = mk_1D_array(m, false);

#pragma omp parallel for num_threads(n_threads) schedule(static)
    for (size_t i = 0; i < m; i++) {
        diag[i] = 2.0 * (1.0 - cos((i+1) * PI / n));
    }

    /*
     * Allocate the matrices b and bt which will be used for storing value of
     * G, \tilde G^T, \tilde U^T, U as described in Chapter 9. page 101.
     */
    real **b = mk_2D_array(m, m, false);
    real **bt = mk_2D_array(m, m, false);

    /*
     * This vector will holds coefficients of the Discrete Sine Transform (DST)
     * but also of the Fast Fourier Transform used in the FORTRAN code.
     * The storage size is set to nn = 4 * n, look at Chapter 9. pages 98-100:
     * - Fourier coefficients are complex so storage is used for the real part
     *   and the imaginary part.
     * - Fourier coefficients are defined for j = [[ - (n-1), + (n-1) ]] while
     *   DST coefficients are defined for j [[ 0, n-1 ]].
     * As explained in the Lecture notes coefficients for positive j are stored
     * first.
     * The array is allocated once and passed as arguments to avoid doings
     * reallocations at each function call.
     */
    int nn = 4 * n;
    real *z[n_threads];

    for (int i = 0; i < n_threads; i++)
    		z[i] = mk_1D_array(nn, false);

    /*
     * Initialize the right hand side data for a given rhs function.
     * Note that the right hand-side is set at nodes corresponding to degrees
     * of freedom, so it excludes the boundary (bug fixed by petterjf 2017).
     *
     */
#pragma omp parallel for num_threads(n_threads) schedule(static) collapse(2)
    for (size_t i = from; i < to; i++) {
        for (size_t j = 0; j < m; j++) {
            b[i][j] = h * h * rhs(grid[i+1], grid[j+1]);
        }
    }

    /*
     * Compute \tilde G^T = S^-1 * (S * G)^T (Chapter 9. page 101 step 1)
     * Instead of using two matrix-matrix products the Discrete Sine Transform
     * (DST) is used.
     * The DST code is implemented in FORTRAN in fsf.f and can be called from C.
     * The array zz is used as storage for DST coefficients and internally for
     * FFT coefficients in fst_ and fstinv_.
     * In functions fst_ and fst_inv_ coefficients are written back to the input
     * array (first argument) so that the initial values are overwritten.
     */

#pragma omp parallel for num_threads(n_threads) schedule(static)
    for (size_t i = from; i < to; i++) {
        fst_(b[i], &n, z[omp_get_thread_num()], &nn);
    }
    transpose(bt, b, m);

#pragma omp parallel for num_threads(n_threads) schedule(static)
    for (size_t i = from; i < to; i++) {
      fstinv_(bt[i], &n, z[omp_get_thread_num()], &nn);
    }

    /*
     * Solve Lambda * \tilde U = \tilde G (Chapter 9. page 101 step 2)
     */
#pragma omp parallel for num_threads(n_threads) schedule(static) collapse(2)
    for (size_t i = from; i < to; i++) {
        for (size_t j = 0; j < m; j++) {
            bt[i][j] = bt[i][j] / (diag[i] + diag[j]);
        }
    }

    /*
     * Compute U = S^-1 * (S * Utilde^T) (Chapter 9. page 101 step 3)
     */
#pragma omp parallel for num_threads(n_threads) schedule(static)
    for (size_t i = from; i < to; i++) {
      fst_(bt[i], &n, z[omp_get_thread_num()], &nn);
    }

    transpose(b, bt, m);

 #pragma omp parallel for num_threads(n_threads) schedule(static)
    for (size_t i = from; i < to; i++) {
      fstinv_(b[i], &n, z[omp_get_thread_num()], &nn);
        //fstinv_(b[i], &n, z, &nn);
    }

    /*
     * Compute maximal value of solution for convergence analysis in L_\infty
     * norm.
     */
    double u_max = 0.0;

// TODO Doesn't these two do the same? One with omp, other with MPI? Check this
#pragma omp parallel for num_threads(n_threads) reduction(max: u_max)
    for (size_t i = from; i < to; i++) {
        for (size_t j = 0; j < m; j++) {
            u_max = u_max > b[i][j] ? u_max : b[i][j];
        }
    }
    double u_max_after_reduce = 0.0;
    MPI_Allreduce(&u_max, &u_max_after_reduce, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if(mpi_rank == 0){
        printf("u_max = %e\n", u_max_after_reduce);
    }

    //Calculating error, from appendix B
    double maxerr = 0.0;

    for (int i = from; i < to; i++)
    {
    	for (int j = 0; j < m; j++)
    	{
    		real error = fabs(answer(grid[i+1], grid[j+1]) - b[i][j]);
            maxerr = maxerr > error ? maxerr : error;
    	}
    }
    double error_u_max = 0.0;
    MPI_Allreduce (&maxerr, &error_u_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if(mpi_rank == 0){
      printf ("Largest error encountered: %e\n", error_u_max);
    }

    // Completed time
    if(mpi_rank == 0){
        double time_final = MPI_Wtime() - time1;
        printf("Time: %f\n", time_final);
    }

    free(sendcounts);
    free(recvcounts);
    free(rdispls);
    free(sdispls);
    free(from_proc);
    free(to_proc);
    MPI_Finalize ();
    return 0;

}

/*
 * This function is used for initializing the right-hand side of the equation.
 * Other functions can be defined to swtich between problem definitions.
 */

real rhs(real x, real y) {
    //return 2 * (y - y*y + x - x*x);
    return 5.0 * PI * PI * sin (PI * x) * sin (2.0 * PI * y);
}

// the exact solution, provided in appendix B
real answer(real x, real y) {
    return sin(PI * x) * sin(2 * PI * y);
}
/*
 * Write the transpose of b a matrix of R^(m*m) in bt.
 * In parallel the function MPI_Alltoallv is used to map directly the entries
 * stored in the array to the block structure, using displacement arrays.
 */
// Denne er viktig, her er alt i loopen i samme prosess, vær nøye her. Må flippes, transponereres
void transpose(real **bt, real **b, size_t m)
{
    MPI_Alltoallv( b[0], sendcounts, sdispls, MPI_DOUBLE, bt[0], recvcounts, rdispls, m_columntype, MPI_COMM_WORLD );
}



/*
 * The allocation of a vectore of size n is done with just allocating an array.
 * The only thing to notice here is the use of calloc to zero the array.
 */

real *mk_1D_array(size_t n, bool zero)
{
    if (zero) {
        return (real *)calloc(n, sizeof(real));
    }
    return (real *)malloc(n * sizeof(real));
}

/*
 * The allocation of the two-dimensional array used for storing matrices is done
 * in the following way for a matrix in R^(n1*n2):
 * 1. an array of pointers is allocated, one pointer for each row,
 * 2. a 'flat' array of size n1*n2 is allocated to ensure that the memory space
 *   is contigusous,
 * 3. pointers are set for each row to the address of first element.
 */

real **mk_2D_array(size_t n1, size_t n2, bool zero)
{
    // 1
    real **ret = (real **)malloc(n1 * sizeof(real *));

    // 2
    if (zero) {
        ret[0] = (real *)calloc(n1 * n2, sizeof(real));
    }
    else {
        ret[0] = (real *)malloc(n1 * n2 * sizeof(real));
    }

    // 3
    for (size_t i = 1; i < n1; i++) {
        ret[i] = ret[i-1] + n2;
    }
    return ret;
}
