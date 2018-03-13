#include <math.h>
#include "zeta0.h"
#include <stdio.h>

// Zeta function
double zeta_function(int n)
{
  double sum = 0;

  for (double i = 1; i <= n; i++){
    sum += 1/(i*i);
  }

  double the_pi = sqrt(sum*6);

  return the_pi;
}
