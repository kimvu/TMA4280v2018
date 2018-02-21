#include <stdio.h>
#include <assert.h>
#include "zeta0.h"
#include <math.h>
#define M_PI 3.14159265358979323846

int verification_zeta0()
{
  double errors[24];
    for (int i = 1; i <= 24; i++){
      errors[i-1] = (fabs(M_PI - zeta_function(pow(2, i))));
    }

    FILE *f = fopen("verification_results.txt", "w");

    for (int i = 1; i <= 24; i++){
      printf("%f\n", errors[i-1]);
      fprintf(f, "%f\n", errors[i-1]);
    }

}

int main(int argc, char *argv[])
{
    int ret = 0;
    ret |= verification_zeta0();
    return ret;
}
