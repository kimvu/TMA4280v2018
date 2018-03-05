#include <stdio.h>
#include <assert.h>
#include "zeta0.h"
#include <math.h>

int test_zeta0()
{
    double results = zeta_function(3);
    printf("Pi is approximately: %f\n", results);
    double expected_value = 2.858;

    double difference = results - expected_value;
    if (fabs(difference) > 0.001){
      perror("The function is not returning the expected value");
      return (-1);
    }
    return 0;
}

int main(int argc, char *argv[])
{
    int ret = 0;
    ret = test_zeta0();
    return ret;
}
