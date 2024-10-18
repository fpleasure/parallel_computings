#include <iostream>
#include <ctime>
#include <cmath>

#include <vector>

const int n = 1024;

void createMatrix(double**& A, double**& B, double**& C) {
	A = new double*[n];
	B = new double* [n];
	C = new double* [n];

	for (int i = 0; i < n; ++i) {
		A[i] = new double[n];
		B[i] = new double[n];
		C[i] = new double[n];
	}
}

void deleteMatrix(double** A, double** B, double** C) {
	for (int i = 0; i < n; ++i) {
		delete[] A[i];
		delete[] B[i];
		delete[] C[i];
	}

	delete[] A;
	delete[] B;
	delete[] C;
}

void initMatrix(double **A, double **B, double **C) {
	for (int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j){
			A[i][j] = sin(i + j);
			B[i][j] = cos(i - j);
			C[i][j] = 0.0;
	}
}

void clearC(double** C) {
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			C[i][j] = 0.0;
}

void multiplyIJK(double **A, double **B, double **C)
{
	clock_t t1, t2;
	t1 = clock();

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			for (int k = 0; k < n; ++k)
				C[i][j] += A[i][k] * B[k][j];

	t2 = clock();
	printf("Time for n=%d (IJK): %f, C[n/2][n/2] = %f\n", n, (double)(t2 - t1) / CLOCKS_PER_SEC, C[n / 2][n / 2]);
}

void multiplyJKI(double** A, double** B, double** C)
{
	clock_t t1, t2;
	t1 = clock();

	for (int j = 0; j < n; ++j)
		for (int k = 0; k < n; ++k)
			for (int i = 0; i < n; ++i)
				C[i][j] += A[i][k] * B[k][j];

	t2 = clock();
	printf("Time for n=%d (JKI): %f, C[n/2][n/2] = %f\n", n, (double)(t2 - t1) / CLOCKS_PER_SEC, C[n / 2][n / 2]);
}

void multiplyKIJ(double** A, double** B, double** C)
{
	clock_t t1, t2;
	t1 = clock();

	for (int k = 0; k < n; ++k)
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				C[i][j] += A[i][k] * B[k][j];

	t2 = clock();
	printf("Time for n=%d (KIJ): %f, C[n/2][n/2] = %f\n", n, (double)(t2 - t1) / CLOCKS_PER_SEC, C[n / 2][n / 2]);
}

void multiplyJIK(double** A, double** B, double** C)
{
	clock_t t1, t2;
	t1 = clock();

	for (int j = 0; j < n; ++j)
		for (int i = 0; i < n; ++i)		
			for (int k = 0; k < n; ++k)
				C[i][j] += A[i][k] * B[k][j];

	t2 = clock();
	printf("Time for n=%d (JIK): %f, C[n/2][n/2] = %f\n", n, (double)(t2 - t1) / CLOCKS_PER_SEC, C[n / 2][n / 2]);
}

void multiplyIKJ(double** A, double** B, double** C)
{
	clock_t t1, t2;
	t1 = clock();

	for (int i = 0; i < n; ++i)
		for (int k = 0; k < n; ++k)
			for (int j = 0; j < n; ++j)			
				C[i][j] += A[i][k] * B[k][j];

	t2 = clock();
	printf("Time for n=%d (IKJ): %f, C[n/2][n/2] = %f\n", n, (double)(t2 - t1) / CLOCKS_PER_SEC, C[n / 2][n / 2]);
}

void multiplyKJI(double** A, double** B, double** C)
{
	clock_t t1, t2;
	t1 = clock();

	for (int k = 0; k < n; ++k)
		for (int j = 0; j < n; ++j)
			for (int i = 0; i < n; ++i)			
				C[i][j] += A[i][k] * B[k][j];

	t2 = clock();
	printf("Time for n=%d (KJI): %f, C[n/2][n/2] = %f\n", n, (double)(t2 - t1) / CLOCKS_PER_SEC, C[n / 2][n / 2]);
}

int main()
{
	double **A, **B, **C;
	
	createMatrix(A, B, C);
	initMatrix(A, B, C);

	multiplyIJK(A, B, C);

	clearC(C);
	multiplyJKI(A, B, C);
	
	clearC(C);
	multiplyKIJ(A, B, C);

	clearC(C);
	multiplyIKJ(A, B, C);

	clearC(C);
	multiplyKJI(A, B, C);

	clearC(C);
	multiplyJIK(A, B, C);
	
	deleteMatrix(A, B, C);

	return EXIT_SUCCESS;
}
