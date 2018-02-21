#include <stdio.h>
#include <assert.h>
#include "mach0.h"
#include <math.h>

int test_mach0()
{
    double results = mach_function(3);
    printf("Pi is approximately: %f\n", results);
    double expected_value = 3.060;

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
    ret = test_mach0();
    return ret;
}
