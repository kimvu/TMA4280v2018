#include "mach2.h"
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>

// The arctan function
double arctan(int n_threads,int n, double x)
{
  double sum = 0;

  #pragma omp parallell num_threads(n_threads) reduction(+:sum)
  for (double i = 1; i <= n; i++){
    sum += pow(-1, i-1)*(pow(x, 2*i-1)/(2*i-1));
  }
  return sum;
}

// The mach function
double mach2_function(int n_threads, int n){
      for (int i=1; i<=n; i++) {
        double first_term = 16*arctan(n_threads, n, (double)0.2);
        double second_term = 4*arctan(n_threads, n, (double)1/239);

        return first_term - second_term;
      }
  }

// Main funciton, tests the code
int main(int argc, char *argv[]){
     int n = atoi(argv[1]);
     int n_threads = atoi(argv[2]);

     double res = mach2_function(n_threads, n);

     printf("Result: %f\n", res);
     exit ( EXIT_SUCCESS );
}
