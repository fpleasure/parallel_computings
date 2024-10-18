#pragma once
#include <iostream>
#include <vector>
#include <chrono>

/*
TODO:

1) 2-й вариант LU decomposition
2) Блочный вариант

*/

std::vector<double> make_LU_decomposition(const std::vector<double> & M, int n) {
	// LU decomposition for matrix A
	std::vector<double> A = M;
	printf("[INFO] Make LU decomposition\n");
	auto start_time = std::chrono::high_resolution_clock::now();
	
	for (int k = 0; k < n; ++k) {
		for (int i = k + 1; i < n; ++i) {
			A[i * n + k] /= A[k * n + k];
			for(int j = k + 1; j < n; ++j) {
				A[i * n + j] -= A[i * n + k] * A[k * n + j];
			}
		}
	}

	auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
	
	printf("[INFO] Success LU decomposition with time = %f seconds\n", elapsed.count());
	return A;
}

/*
std::vector<double> get_block_for_block_LU_decomposition(std::vector<double> & M, int n) {

}

std::vector<double> make_block_LU_decomposition(std::vector<double> & M, int n) {

}
*/

void print_matrix(const std::vector<double> & A, int n) {
	// Print matrix, saving as vector
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			std::cout << A[n * i + j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

std::vector<double> cross_matrix(const std::vector<double> & A, const std::vector<double> & B, int n) {
	// Cross matrix, saving as vector
	std::vector<double> C = std::vector<double>(n * n, 0);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			for (int k = 0; k < n; ++k)
				C[i * n + j] += A[i * n + k] * B[k * n + j];
	return C;
}

std::vector<double> get_L_from_LU_decomposition(const std::vector<double> & LU, int n) {
	// Get L from LU decomposition
	std::vector<double> L = std::vector<double>(n * n, 1);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (i > j) {
				L[n * i + j] = LU[n * i + j];
			}
			else if (i == j) {
				L[n * i + j] = 1;
			}
			else {
				L[n * i + j] = 0;
			}
		}
	}
	return L;
}

std::vector<double> get_U_from_LU_decomposition(const std::vector<double> & LU, int n) {
	// Get U from LU decomposition
	std::vector<double> U = std::vector<double>(n * n, 1);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (i <= j) {
				U[n * i + j] = LU[n * i + j];
			}
			else {
				U[n * i + j] = 0;
			}
		}
	}
	return U;
}

bool test_LU_decomposition(const std::vector<double> & LU, const std::vector<double> & A, int n) {
	// Test for LU decomposition
	printf("[INFO] Matrix L from LU decomposition:\n");
	std::vector<double> L = get_L_from_LU_decomposition(LU, n);
	print_matrix(L, n);
	printf("[INFO] Matrix U from LU decomposition:\n");
	std::vector<double> U = get_U_from_LU_decomposition(LU, n);
	print_matrix(U, n);
	printf("[INFO] Matrix L*U:\n");
	std::vector<double> C = cross_matrix(L, U, n);
	print_matrix(C, n);
	printf("[INFO] First matrix:\n");
	print_matrix(A, n);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (C[i * n +j] != A[i * n +j]) {
				printf("[ERROR] Matrix L*U not equals A!\n");
				return false;
			}
		}
	}
	printf("[INFO] Test LU decomposition passed\n");
	return true;
	
}
