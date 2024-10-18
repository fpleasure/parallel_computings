// VectorSum.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include <iostream>

#include <vector>
#include <cmath>

#include "omp.h"

double f(double x) {
    return pow(pow(sqrt(sqrt(fabs(sin(log(x))))), 0.98), 0.76);
}

int main()
{
    const int np = omp_get_max_threads();
    printf("Number of threads: %d\n", np);

    const int n = 10000000;
    std::vector<double> a(n);
    //double* a = new double[n];
    for (int i = 0; i < n; ++i)
        a[i] = 1.0 + (double)rand() / RAND_MAX;

    double t1, t2;

    //последовательный режим
    double sum(0.0);

    t1 = omp_get_wtime();
    for (int i = 0; i < n; ++i)
        sum += f(a[i]);
    t2 = omp_get_wtime();

    printf("Sequential mode: sum = %f, time = %f\n", sum, t2 - t1);

    //параллельный режим 1 (пример гонки данных)
    sum = 0.0;

    t1 = omp_get_wtime();
#pragma omp parallel
    {
        const int id = omp_get_thread_num();
        //printf("Thread %d\n", id);
        
        for (int i = 0; i < n; ++i)
            sum += f(a[i]);
    }
    t2 = omp_get_wtime();

    printf("Parallel mode 1 (incorrect): sum = %f, time = %f\n", sum, t2 - t1);

    //параллельный режим 2 (ручное разделение работы)
    sum = 0.0;
    std::vector<double> partialSums(np);
    //double* partialSums = new double[np];


    t1 = omp_get_wtime();
#pragma omp parallel
    {
        const int id = omp_get_thread_num();
        //printf("Thread %d\n", id);

        partialSums[id] = 0.0;

        for (int i = id * n / np; i < (id + 1) * n / np; ++i)
            partialSums[id] += f(a[i]);
    }
    t2 = omp_get_wtime();

    for (int i = 0; i < np; ++i)
        sum += partialSums[i];

    printf("Parallel mode 2 (manual workload division): sum = %f, time = %f\n", sum, t2 - t1);

    //параллельный режим 3 (ручное разделение работы)
    sum = 0.0;

    t1 = omp_get_wtime();
#pragma omp parallel
    {
        const int id = omp_get_thread_num();
        //printf("Thread %d\n", id);

        partialSums[id] = 0.0;

        for (int i = id; i < n; i += np)
            partialSums[id] += f(a[i]);
    }
    t2 = omp_get_wtime();

    for (int i = 0; i < np; ++i)
        sum += partialSums[i];

    printf("Parallel mode 3 (manual workload division): sum = %f, time = %f\n", sum, t2 - t1);

    //параллельный режим 4 (критическая секция, неэффективно)
    sum = 0.0;

    t1 = omp_get_wtime();
#pragma omp parallel
    {
        const int id = omp_get_thread_num();
        //printf("Thread %d\n", id);

        for (int i = id; i < n; i += np)
#pragma omp critical            
            sum += f(a[i]);
    }
    t2 = omp_get_wtime();

    printf("Parallel mode 4 (critical section - ineffective): sum = %f, time = %f\n", sum, t2 - t1);

    //параллельный режим 5 (критическая секция, эффективно)
    sum = 0.0;

    t1 = omp_get_wtime();
#pragma omp parallel
    {
        const int id = omp_get_thread_num();
        //printf("Thread %d\n", id);
        double localSum(0.0);

        for (int i = id; i < n; i += np)
            localSum += f(a[i]);

#pragma omp critical
            sum += localSum;
    }
    t2 = omp_get_wtime();

    printf("Parallel mode 5 (critical section - effective): sum = %f, time = %f\n", sum, t2 - t1);

    //параллельный режим 6 (atomic-операции, неэффективно)
    sum = 0.0;

    t1 = omp_get_wtime();
#pragma omp parallel
    {
        const int id = omp_get_thread_num();
        //printf("Thread %d\n", id);

        for (int i = id; i < n; i += np)
#pragma omp atomic            
            sum += f(a[i]);
    }
    t2 = omp_get_wtime();

    printf("Parallel mode 6 (atomic operations - ineffective): sum = %f, time = %f\n", sum, t2 - t1);

    //параллельный режим 7 (atomic-операции, эффективно)
    sum = 0.0;

    t1 = omp_get_wtime();
#pragma omp parallel
    {
        const int id = omp_get_thread_num();
        //printf("Thread %d\n", id);
        double localSum(0.0);

        for (int i = id; i < n; i += np)
            localSum += f(a[i]);

#pragma omp atomic
        sum += localSum;
    }
    t2 = omp_get_wtime();

    printf("Parallel mode 7 (atomic section - effective): sum = %f, time = %f\n", sum, t2 - t1);

    //параллельный режим 8 (замки)
    sum = 0.0;
    omp_lock_t lock;
    omp_init_lock(&lock);

    t1 = omp_get_wtime();
#pragma omp parallel shared(lock,a,sum) default(none)
    {
        const int id = omp_get_thread_num();
        //printf("Thread %d\n", id);
        double localSum(0.0);

        for (int i = id; i < n; i += np)
            localSum += f(a[i]);

        omp_set_lock(&lock);
        sum += localSum;
        omp_unset_lock(&lock);
    }
    t2 = omp_get_wtime();

    printf("Parallel mode 8 (locks): sum = %f, time = %f\n", sum, t2 - t1);

    omp_destroy_lock(&lock);

    //параллельный режим 9 (редукция)
    sum = 0.0;

    t1 = omp_get_wtime();
#pragma omp parallel reduction(+: sum)
    {
        const int id = omp_get_thread_num();
        //printf("Thread %d\n", id);

        for (int i = id; i < n; i += np)
            sum += f(a[i]);
    }
    t2 = omp_get_wtime();

    printf("Parallel mode 9 (reduction): sum = %f, time = %f\n", sum, t2 - t1);

    //параллельный режим 10 (редукция)
    sum = 0.0;

    t1 = omp_get_wtime();
#pragma omp parallel for reduction(+: sum)
    for (int i = 0; i < n; ++i)
        sum += f(a[i]);
    t2 = omp_get_wtime();

    printf("Parallel mode 10 (reduction): sum = %f, time = %f\n", sum, t2 - t1);


    //delete[] a;
    //delete[] partialSums;

    return EXIT_SUCCESS;
}

