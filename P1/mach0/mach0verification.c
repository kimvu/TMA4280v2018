#include <stdio.h>
#include <assert.h>
#include "mach0.h"
#include <math.h>
#define M_PI 3.14159265358979323846

int verification_mach0()
{
  double errors[24];
    for (int i = 1; i <= 24; i++){
      errors[i-1] = (fabs(M_PI - mach_function(pow(2, i))));
    }

    FILE *f = fopen("verification_results.txt", "w");

    for (int i = 1; i <= 24; i++){
      printf("%e\n", errors[i-1]);
      fprintf(f, "%e\n", errors[i-1]);
    }
    return 0;
}

int main(int argc, char *argv[])
{
    int ret = 0;
    ret |= verification_mach0();
    return ret;
}
