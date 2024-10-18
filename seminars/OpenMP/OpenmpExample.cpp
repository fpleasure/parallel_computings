// OpenmpExample.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <omp.h>

int main()
{
    //omp_set_num_threads(4);
    //omp_set_nested(1);             //MS VC++ (Intel считается устаревшей)
    //omp_set_max_active_levels(2);//только Intel

    int maxThreads = omp_get_max_threads();

    int numThreads = omp_get_num_threads();
    printf("Maximum/actual number of threads = %d/%d\n", maxThreads, numThreads);

    int i, j;
    std::cin >> i >> j;

    int commonVariable = 100;
    int variable2 = 250;

    printf("Address of variable2 = %p, value = %d\n", &variable2, variable2);

//#pragma omp parallel num_threads(i) shared(commonVariable,std::cout) firstprivate(variable2) default(none)
#pragma omp parallel num_threads(i) firstprivate(variable2)
    {
        int id = omp_get_thread_num();
        printf("1.1 Thread %d (address of id = %p) out of %d\n", id, &id, omp_get_num_threads());
        
//        std::cout << id;

        printf("Address of commonVariable = %p, value = %d (thread %d)\n",
            &commonVariable, commonVariable, id);

        variable2 += id;

        printf("Address of variable2 = %p, value = %d (thread %d)\n",
            &variable2, variable2, id);

        //вложенный параллелизм
        /*
#pragma omp parallel num_threads(2)
        {
            int id = omp_get_thread_num();
            printf("1.2 Thread %d out of %d\n", id, omp_get_num_threads());
        }*/

        /*std::cout << "Hello ";
        std::cout << "World!";
        std::cout << std::endl;*/
    }

    printf("Address of variable2 = %p, value = %d\n", &variable2, variable2);

#pragma omp parallel num_threads(j)
//#pragma omp parallel if(j>=4)
    {
        int id = omp_get_thread_num();
        printf("2. Thread %d out of %d\n", id, omp_get_num_threads());
    }


    return EXIT_SUCCESS;
}
