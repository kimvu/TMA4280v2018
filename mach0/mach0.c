#include "mach0.h"

#include <stdio.h>
#include <math.h>

void arctan(int n, double x)
{
  double sum = 0;

  for (double i = 1; i <= n; i++){
    sum += pow(-1, i-1)*(pow(x, 2*i-1)/(2*1-1));
  }

  return sum;
}

void mach_function(int n){
  double first_term = 16*arctan(n, (double)0.2);
  double second_term = 4*arctan(n, (double)1/239);

  double the_pi = first_term - second_term;

  return the_pi;
}
