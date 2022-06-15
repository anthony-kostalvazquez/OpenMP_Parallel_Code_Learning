#include <stdio.h>
#include "omp.h"

#define NUM_THREADS 12
#define PAD 8
static long num_steps = 1000000000;

double non_para(long);

double para(long);

double para_padded(long);

double para_synch(long);

int main()
{
    double start, end; 

    //runs non parallel code for the integration of f(x) = 4/(1+x^2) from 0 to 1
    start = omp_get_wtime();
    printf("The non parallel code got the solution %f \n", non_para(num_steps));
    end = omp_get_wtime();
    printf("Time: %f \n", (end - start));

    //runs parallel code for the integration of f(x) = 4/(1+x^2) from 0 to 1
    start = omp_get_wtime();
    printf("The spmd parallel code got the solution %f \n", para(num_steps));
    end = omp_get_wtime();
    printf("Time: %f \n", (end - start));

    //runs parallel code for the integration of f(x) = 4/(1+x^2) from 0 to 1
    start = omp_get_wtime();
    printf("The spmd parallel code with the padding got the solution %f \n", para_padded(num_steps));
    end = omp_get_wtime();
    printf("Time: %f \n", (end - start));

    //runs parallel code for the integration of f(x) = 4/(1+x^2) from 0 to 1
    start = omp_get_wtime();
    printf("The sychronized parallel code got the solution %f \n", para_padded(num_steps));
    end = omp_get_wtime();
    printf("Time: %f \n", (end - start));

}

//integrates f(x) = 4/(1+x^2) from 0 to 1 and returns pi
double non_para(long num_steps)
{
    double step, x, pi, sum = 0.0;
    step = 1.0 / num_steps;

    //this for loop will calculate the integral of f(x) = 4/(1+x^2) from 0 to 1
    for(int i = 0; i < num_steps; i++)
    {
        x = (i + 0.5) * step;
        sum = sum + 4.0/(1.0 + x*x);
    }
    pi = step * sum;

    return pi;
}

double para(long num_steps)
{
    double step, pi, sum[NUM_THREADS];
    int i, nthreads;

    step = 1.0 / num_steps;
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel
    {
        int i, id, nthrds;
        double x;
        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();
        
        //this will ensure that we have the same number of threads... we ask the system for x amount of threads,  however we dont always get it
        if (id == 0)
        {
            nthreads = nthrds;
        }

        for(i = id, sum[id] = 0.0 ; i < num_steps; i = i + nthrds)
        {
            x = (i + 0.5) * step;
            sum[id] +=  + 4.0/(1.0 + x*x);
        }

        //suming the array of each threads sum
        for(i = 0, pi = 0.0; i < nthreads; i++)
        {
            pi += sum[i] * step;
        }
    }



    return pi;
}

//padding requires alot of knowledge of the hardware
double para_padded(long num_steps)
{
    double step, pi, sum[NUM_THREADS][PAD];
    int i, nthreads;

    step = 1.0 / num_steps;
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel
    {
        int i, id, nthrds;
        double x;
        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();
        
        //this will ensure that we have the same number of threads... we ask the system for x amount of threads,  however we dont always get it
        if (id == 0)
        {
            nthreads = nthrds;
        }

        for(i = id, sum[id][0] = 0.0 ; i < num_steps; i = i + nthrds)
        {
            x = (i + 0.5) * step;
            sum[id][0] +=  + 4.0/(1.0 + x*x);
        }

        //suming the array of each threads sum
        for(i = 0, pi = 0.0; i < nthreads; i++)
        {
            pi += sum[i][0] * step;
        }
    }



    return pi;
}

//utilizes a critical statement to add at the end. this will be basicly as fast as the padded spmd 
double para_synch(long num_steps)
{
    double step, pi;
    int i, nthreads;

    step = 1.0 / num_steps;

    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel
    {
        int i, id, nthrds;
        double x, sum;
        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();
        
        //this will ensure that we have the same number of threads... we ask the system for x amount of threads,  however we dont always get it
        if (id == 0)
        {
            nthreads = nthrds;
        }

        for(i = id, sum = 0.0 ; i < num_steps; i = i + nthrds)
        {
            x = (i + 0.5) * step;
            sum +=  + 4.0/(1.0 + x*x);
        }

        #pragma omp critical
        {
            pi += sum * step;
        }
    }



    return pi;
}