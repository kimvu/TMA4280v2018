#include <stdio.h>
#include <assert.h>
#include "zeta0.h"

void test_zeta0()
{
    double results = zeta_function(3);
    printf("Pi is approximately: " + "%d\n", results);
}


int main(int argc, char *argv[])
{
    test_zeta0();
    return 0
}
