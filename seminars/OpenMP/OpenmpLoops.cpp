// OpenmpLoops.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <fstream>
#include <vector>
#include "omp.h"

const int n = 100000;

double f(double x) {
    return pow(pow(sqrt(sqrt(fabs(sin(log(x))))), 0.98), 0.76);
}

void taskFunc(int taskNum) {
    int count = (rand() % 100) + 1;
    double sum = 0.0;
    for (int j = 0; j < count; ++j) {
        double x = 1.0 + (double)rand() / RAND_MAX;
        sum += f(x);
    }
        
    printf("Task %d, thread %d, sum = %f\n", taskNum, omp_get_thread_num(), sum);
}

double x;
#pragma omp threadprivate(x)

int main()
{
    omp_set_num_threads(6);

    const int np = omp_get_max_threads();
    printf("Number of threads: %d\n", np);

    srand(time(NULL));


    std::vector<double> a(n);
    for (int i = 0; i < n; ++i)
        a[i] = 1.0 + (double)rand() / RAND_MAX;

    double time = -omp_get_wtime();

    /*int i = 1;
    for ( ;  ; ) {

        i += 2;

        if (i >= n)
            break;
    }*/

    std::vector<unsigned int> elementThread(n);

    int q;

    //volatile double res = 0.0;
#pragma omp parallel
    {
        volatile double res = 0.0;
        int totalCount = 0;

        const int id = omp_get_thread_num();

//#pragma omp for schedule(guided, 100)
#pragma omp for lastprivate(q)
        for (int i = 0; i < n; ++i) {
            int count = (rand() % 1000) + 1;
            for (int j = 0; j < count; ++j) {
                double x = 1.0 + (double)rand() / RAND_MAX;
                res += f(x);
            }

            totalCount += count;

            elementThread[i] = id + 1;

            q = i;
        }

        printf("Thread %d, sum = %f (%d numbers)\n", id, res, totalCount);
    }

    time += omp_get_wtime();
    //printf("Time = %f, res = %f\n", time, res);
    printf("Time = %f, q = %d\n", time, q);

    std::ofstream file("iterations.dat");
    if (file.is_open()) {
        for (int i = 0; i < n; ++i)
            file << i + 1 << " " << elementThread[i] << std::endl;

        file.close();
    }

#pragma omp parallel for ordered
    for (int i = 0; i < 10; ++i) {
#pragma omp ordered
        printf("Iteration %d, thread %d\n", i, omp_get_thread_num());
    }

#pragma omp parallel
    {
        const int id = omp_get_thread_num();
        double sum = 0.0;

#pragma omp for nowait//collapse(2) //не работает с MS VC++
        for (int i = 0; i < 10; ++i)
            for (int j = 0; j < 100000; ++j)
                sum += f(a[j]);
        //sum += f(a[i * 10000 + j]);

/*#pragma omp for
        for (int k = 0; k < 10 * 10000; ++k) {
            int i = k / 10000;
            int j = k % 10000;
        }*/

        printf("1. Thread %d, sum = %f\n", id, sum);

        sum = 0.0;

#pragma omp for 
        for (int i = 0; i < 10; ++i)
            for (int j = 0; j < 1000; ++j)
                sum += f(a[i * 10000 + j]);

        printf("2. Thread %d, sum = %f\n", id, sum);

#pragma omp sections
        {
#pragma omp section
            {
                volatile double sum = 0.0;
                for (int i = 0; i < n; ++i)
                    sum += f(a[i]);
                printf("Section A - Thread %d, sum = %f\n", omp_get_thread_num(), sum);
            }

#pragma omp section
            {
                volatile double sum = 0.0;
                for (int i = 0; i < n; ++i)
                    sum += f(a[i]);
                printf("Section B - Thread %d, sum = %f\n", omp_get_thread_num(), sum);
            }

#pragma omp section
            {
                volatile double sum = 0.0;
                for (int i = 0; i < n; ++i)
                    sum += f(a[i]);
                printf("Section C - Thread %d, sum = %f\n", omp_get_thread_num(), sum);
            }

#pragma omp section
            {
                volatile double sum = 0.0;
                for (int i = 0; i < n; ++i)
                    sum += f(a[i]);
                printf("Section D - Thread %d, sum = %f\n", omp_get_thread_num(), sum);
            }
        }
    }

#pragma omp parallel
    {
#pragma omp master
        for (int i = 0; i < 10; ++i)
#pragma omp task untied
            taskFunc(i);
    }

    x = 10.0;
    printf("0. x = %f\n", x);

#pragma omp parallel copyin(x)
    {
        const int id = omp_get_thread_num();
        x += 100 * (id + 1);
        printf("1. x = %f, thread = %d\n", x, id);
    }

    printf("2. x = %f\n", x);

#pragma omp parallel
    {
        const int id = omp_get_thread_num();
        //x += 100 * (id + 1);
        printf("3. x = %f, thread = %d\n", x, id);
    }

    return EXIT_SUCCESS;
}
