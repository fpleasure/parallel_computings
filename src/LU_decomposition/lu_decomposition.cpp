#pragma once

#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <omp.h>
#include "matrix.cpp"

/**
 * LU-разложение матрицы M размера m*n, m >= n. Результат записывается в M.
 */
double make_LU_decomposition(std::vector<double>& A, int m, int n) {
    printf("[INFO] Make LU decomposition\n");
    auto start_time = std::chrono::high_resolution_clock::now();

    for (int k = 0; k < std::min(m - 1, n); ++k) {
        for (int i = k + 1; i < m; ++i) {
            A[i * n + k] /= A[k * n + k];

            if (k < n) {
                for (int j = k + 1; j < n; ++j) {
                    A[i * n + j] -= A[i * n + k] * A[k * n + j];
                }
            }
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    printf("[INFO] Success LU decomposition with time = %f seconds\n", elapsed.count());

    return elapsed.count();
}

/**
 * Параллельное LU-разложение матрицы M размера m*n, m >= n. Результат записывается в M.
 */
double make_parallel_LU_decomposition(std::vector<double>& A, int m, int n) {
    printf("[INFO] Make parallel LU decomposition\n");
    auto start_time = std::chrono::high_resolution_clock::now();

    for (int k = 0; k < std::min(m - 1, n); ++k) {
        // Первая фаза: деление столбца A[i * n + k] на главный элемент A[k * n + k]
#pragma omp parallel for
        for (int i = k + 1; i < m; ++i) {
            A[i * n + k] /= A[k * n + k];
        }

        // Вторая фаза: обновление подматрицы
#pragma omp parallel for
        for (int i = k + 1; i < m; ++i) {
            for (int j = k + 1; j < n; ++j) {
                A[i * n + j] -= A[i * n + k] * A[k * n + j];
            }
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    printf("[INFO] Success parallel LU decomposition with time = %f seconds\n", elapsed.count());

    return elapsed.count();
}

/**
 * LU-разложение подматрицы размера sm*sn, sm >= sn, исходной матрицы M размера m*n,
 * левый верхний элемент которой имеет индексы (offsetI, offsetJ).
 */
double make_sub_LU_decomposition(std::vector<double>& A, int m, int n, 
                                 int sm, int sn, int offsetI, int offsetJ) {
    printf("[INFO] Make sub LU decomposition\n");

    auto start_time = std::chrono::high_resolution_clock::now();

    for (int k = 0; k < std::min(sm - 1, sn); ++k) {
        for (int i = offsetI + k + 1; i < offsetI + sm; ++i) {
            A[i * n + k] /= A[k * n + k];
            if (k < sn) {
                for (int j = offsetJ + k + 1; j < offsetJ + sn; ++j) {
                    A[i * n + j] -= A[i * n + k] * A[k * n + j];
                }
            }
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    printf("[INFO] Success sub LU decomposition with time = %f seconds\n", elapsed.count());

    return elapsed.count();
}

/**
 * Блочное LU-разложение матрицы размера n*n с блоком размера b.
 */
double make_block_LU_decomposition(std::vector<double>& A, int n, int b) {
    printf("[INFO] Make block LU decomposition\n");
    auto start_time = std::chrono::high_resolution_clock::now();

    // Выделение памяти для блока и полосы под ним
    std::vector<double> B(n * b, 0);
    // Выделение памяти для полосы справа от блока
    int l0 = n - b;
    std::vector<double> C(b * l0, 0);

    for (int k = 0; k < n; k += b) {
        // Копирование блока и полосы под ним в B
        for (int i = 0; i < n - k; ++i) {
            for (int j = 0; j < b; ++j) {
                B[i * b + j] = A[(i + k) * n + j + k];
            }
        }

        // Копирование полосы справа от блока в C
        int l = n - b - k;
        for (int i = 0; i < b; ++i) {
            for (int j = 0; j < l; ++j) {
                C[i * l0 + j] = A[(i + k) * n + j + k + b];
            }
        }

        make_LU_decomposition(B, n - k, b);

        // Вычисление U_23
        for (int i = 0; i < b; ++i) {
            for (int p = 0; p < l; ++p) {
                for (int j = 0; j < i; ++j) {
                    C[i * l0 + p] -= B[i * b + j] * C[j * l0 + p];
                }
            }
        }

        // Обновление оставшихся элементов
        for (int i = 0; i < l; i++) {
            for (int s = 0; s < b; s++) {
                for (int j = 0; j < l; j++) {
                    A[(i + k + b) * n + j + k + b] -= B[(i + b) * b + s] * C[s * l0 + j];
                }
            }
        }

        // Обновление элементов исходной матрицы
        for (int i = 0; i < n - k; i++) {
            for (int j = 0; j < b; j++) {
                A[(i + k) * n + j + k] = B[i * b + j];
            }
        }

        for (int i = 0; i < b; ++i) {
            for (int j = 0; j < l; ++j) {
                A[(i + k) * n + j + b + k] = C[i * l0 + j];
            }
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    printf("[INFO] Success block LU decomposition with time = %f seconds\n", elapsed.count());

    return elapsed.count();
}

/**
 * Параллельное блочное LU-разложение матрицы размера n*n с блоком размера b.
 */
double make_parallel_block_LU_decomposition(std::vector<double>& A, int n, int b) {
    printf("[INFO] Make parallel block LU decomposition\n");
    auto start_time = std::chrono::high_resolution_clock::now();

    // Выделение памяти для блока и полосы под ним
    std::vector<double> B(n * b, 0);
    int l0 = n - b;
    std::vector<double> C(b * l0, 0);

    for (int k = 0; k < n; k += b) {
        // Копирование блока и полосы под ним в B
#pragma omp parallel for collapse(2)
        for (int i = 0; i < n - k; ++i) {
            for (int j = 0; j < b; ++j) {
                B[i * b + j] = A[(i + k) * n + j + k];
            }
        }

        // Копирование полосы справа от блока в C
        int l = n - b - k;
#pragma omp parallel for collapse(2)
        for (int i = 0; i < b; ++i) {
            for (int j = 0; j < l; ++j) {
                C[i * l0 + j] = A[(i + k) * n + j + k + b];
            }
        }

        make_parallel_LU_decomposition(B, n - k, b);

        // Вычисление U_23
#pragma omp parallel for
        for (int p = 0; p < l; ++p) {
            for (int i = 0; i < b; ++i) {
                for (int j = 0; j < i; ++j) {
                    C[i * l0 + p] -= B[i * b + j] * C[j * l0 + p];
                }
            }
        }

        // Обновление оставшихся элементов
#pragma omp parallel for collapse(2)
        for (int i = 0; i < l; i++) {
            for (int j = 0; j < l; j++) {
                for (int s = 0; s < b; s++) {
                    A[(i + k + b) * n + j + k + b] -= B[(i + b) * b + s] * C[s * l0 + j];
                }
            }
        }

        // Обновление элементов исходной матрицы
#pragma omp parallel for collapse(2)
        for (int i = 0; i < n - k; i++) {
            for (int j = 0; j < b; j++) {
                A[(i + k) * n + j + k] = B[i * b + j];
            }
        }

#pragma omp parallel for collapse(2)
        for (int i = 0; i < b; ++i) {
            for (int j = 0; j < l; ++j) {
                A[(i + k) * n + j + b + k] = C[i * l0 + j];
            }
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    printf("[INFO] Success parallel block LU decomposition with time = %f seconds\n", elapsed.count());

    return elapsed.count();
}

/**
 * Вычисление матрицы L из LU разложения.
 */
std::vector<double> get_L_from_LU_decomposition(const std::vector<double>& LU, int n) {
    std::vector<double> L(n * n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i > j) {
                L[n * i + j] = LU[n * i + j];
            } else if (i == j) {
                L[n * i + j] = 1.0;
            }
        }
    }

    return L;
}

/**
 * Вычисление матрицы U из LU разложения.
 */
std::vector<double> get_U_from_LU_decomposition(const std::vector<double>& LU, int n) {
    std::vector<double> U(n * n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i <= j) {
                U[n * i + j] = LU[n * i + j];
            }
        }
    }

    return U;
}

/**
 * Проверка корректности вычислений LU разложения.
 */
bool test_LU_decomposition(const std::vector<double>& LU, const std::vector<double>& A, int n) {
    printf("[INFO] Testing LU decomposition\n");
    std::vector<double> L = get_L_from_LU_decomposition(LU, n);
    std::vector<double> U = get_U_from_LU_decomposition(LU, n);
    std::vector<double> C = cross_matrix(L, U, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (fabs(C[i * n + j] - A[i * n + j]) > 1e-8) {
                printf("[ERROR] The element of L*U differs in modulus from the element of A by more than 1e-8!\n");
                return false;
            }
        }
    }
    
    printf("[INFO] Test LU decomposition passed\n");
    return true;
}
