#pragma once

#include <string>
#include <vector>
#include <fstream>

#include "file.cpp"
#include "matrix.cpp"
#include "lu_decomposition.cpp"

/**
 * Проводит тест с заданным количеством потоков num_threads для матрицы
 * размера n x n, записывая время работы алгоритма в файл с 
 * дескриптором file.
 */
void make_test_LU(int num_threads, int n, std::ofstream& file, 
                  double (*make_method)(std::vector<double>&, int, int), 
                  int block_size = -1) {
    
    printf("[INFO] Make test for LU decomposition: %d threads, matrix n=%d\n", num_threads, n);
    
    omp_set_num_threads(num_threads);
    std::vector<double> B(n * n);
    init_random(B, n);
    
    double time;

    if (block_size == -1) {
        time = make_method(B, n, n);
    } 
    else if (block_size > 0) {
        time = make_method(B, n, block_size);
    }
    else {
        printf("[ERROR] Incorrect block size (=%d) <= 0.\n", block_size);
        return; // Добавляем возврат, чтобы избежать записи времени при ошибке
    }
    
    write_time_to_txt(file, time);
}

/**
 * Проводит несколько тестов make_test_LU, записывая получившуюся 
 * скорость в файл с именем filename; в векторе n_vec задаются
 * размеры матриц.
 */
void make_all_tests_LU(int num_threads, const std::string& filename, 
                        const std::vector<int>& n_vec, 
                        double (*make_method)(std::vector<double>&, int, int), 
                        int block_size = -1) {
    
    printf("[INFO] Make all tests for LU decomposition: %d threads\n", num_threads);
    
    omp_set_num_threads(num_threads);
    std::ofstream file(filename);

    for (const int n : n_vec) {
        make_test_LU(num_threads, n, file, make_method, block_size);
    }
}
