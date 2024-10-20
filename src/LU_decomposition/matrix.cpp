#pragma once

#include <vector>
#include <iostream>
#include <ctime>

/**
 * Выводит матрицу в консоль.
 */
void print_matrix(const std::vector<double>& A, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << A[n * i + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

/**
 * Умножает две матрицы.
 */
std::vector<double> cross_matrix(const std::vector<double>& A, const std::vector<double>& B, int n) {
    std::vector<double> C(n * n, 0.0);
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                C[i * n + j] += A[i * n + k] * B[k * n + j];
            }
        }
    }
    
    return C;
}

/**
 * Инициализирует матрицу случайными значениями.
 */
void init_random(std::vector<double>& A, int n) {
    std::srand(std::time(nullptr));
    
    for (int i = 0; i < n * n; ++i) {
        A[i] = std::rand() % 100 + 0.01; // Генерация случайного числа от 0.01 до 99.99
    }
}