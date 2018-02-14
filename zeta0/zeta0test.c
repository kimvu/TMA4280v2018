#include <stdio.h>
#include <assert.h>
#include "zeta0.h"

void test_zeta0()
{
    double results = zeta_function(3);
    printf("Pi is approximately: " + "%d\n", results);
    double expected_value = 2.858;

    double difference = results - expected_value
    if (fabs(difference) > 0.001){
      perror("The function is not returning the expected value");
      return (-1)
    }
}


int main(int argc, char *argv[])
{
    test_zeta0();
    return 0
}
