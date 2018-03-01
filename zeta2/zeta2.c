#include <math.h>
#include "zeta2.h"
#include <stdio.h>
#include "omp.h"
#include <stdlib.h>


double zeta2_function(int n, int n_threads)
{
  double sum = 0;

  #pragma omp parallell for num_threads(n_threads) reduction(+:sum)
  for (double i = 1; i <= n; i++){
    sum += 1/((double)i*(double)i);
  }

  double the_pi = sqrt(sum*6);

  return the_pi;
}

int main(int argc, char *argv[])
{
   int n = atoi(argv[1]);
   int n_threads = atoi(argv[2]);

   double res = zeta2_function(n, n_threads);

   printf("Result: %f\n", res);
   exit ( EXIT_SUCCESS );
}
