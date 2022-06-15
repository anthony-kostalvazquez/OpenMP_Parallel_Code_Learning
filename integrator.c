#include <stdio.h>
#include "omp.h"

static long num_steps = 100000;
double step; 

int main()
{
    double x, pi, sum = 0.0; 

    step = 1.0 / num_steps;

    //this for loop will calculate the integral of f(x) = 4/(1+x^2) from 0 to 1
    for(int i = 0; i < num_steps; i++)
    {
        x = (i + 0.5) * step;
        sum = sum + 4.0/(1.0 + x*x);
    }
    pi = step * sum;

    printf("%f", pi);
}