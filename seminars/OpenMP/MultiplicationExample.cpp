#include <iostream>
#include <ctime>
#include <cmath>

#include <vector>

#include "mkl.h"

const int n = 8192;

void initMatrix(std::vector<double>& A, std::vector<double>& B, std::vector<double>& C) {
	for (int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j){
			A[i * n + j] = sin(i + j);
			B[i * n + j] = cos(i - j);
			C[i * n + j] = 0.0;
	}
}

void clearC(std::vector<double>& C) {
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			C[i * n + j] = 0.0;
}

void blockMultiply(const double *A, const double *B, double *C, int size, bool measureTime = true)
{
	clock_t t1, t2;
	
	if(measureTime)
		t1 = clock();

	const int bs = 128;
	double a[bs * bs], b[bs * bs], c[bs * bs];

	for(int bi = 0; bi < size; bi += bs)
		for (int bj = 0; bj < size; bj += bs) {
			for (int p = 0; p < bs; ++p)
				for (int q = 0; q < bs; ++q)
					c[p * bs + q] = 0.0;

			for (int bk = 0; bk < size; bk += bs) {
				for (int p = 0; p < bs; ++p)
					for (int q = 0; q < bs; ++q) {
						a[p * bs + q] = A[(bi + p) * size + (bk + q)];
						b[p * bs + q] = B[(bk + p) * size + (bj + q)];
					}

				for (int i = 0; i < bs; ++i)
					for (int k = 0; k < bs; ++k)
						for (int j = 0; j < bs; ++j)
							c[i * bs + j] += a[i * bs + k] * b[k * bs + j];
			}

			for (int p = 0; p < bs; ++p)
				for (int q = 0; q < bs; ++q)
					C[(bi + p) * size + (bj + q)] = c[p * bs + q];
		}

	if (measureTime) {
		t2 = clock();
		printf("Time for block multiplication n=%d (block size = %d): %f, C[n/2][n/2] = %f\n", size, bs, (double)(t2 - t1) / CLOCKS_PER_SEC, C[n / 2 * n + n / 2]);
	}
}

void blockUnrolledMultiply(const std::vector<double>& A, const std::vector<double>& B, std::vector<double>& C)
{
	double t1, t2;
	t1 = clock();
	const int bs = 4;

	double a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33;
	double b00, b01, b02, b03, b10, b11, b12, b13, b20, b21, b22, b23, b30, b31, b32, b33;
	double c00, c01, c02, c03, c10, c11, c12, c13, c20, c21, c22, c23, c30, c31, c32, c33;

	for (int bi = 0; bi < n; bi += bs)
		for (int bj = 0; bj < n; bj += bs)
		{
			c00 = 0.0, c01 = 0.0, c02 = 0.0, c03 = 0.0;
			c10 = 0.0, c11 = 0.0, c12 = 0.0, c13 = 0.0;
			c20 = 0.0, c21 = 0.0, c22 = 0.0, c23 = 0.0;
			c30 = 0.0, c31 = 0.0, c32 = 0.0, c33 = 0.0;

			for (int bk = 0; bk < n; bk += bs)
			{
				a00 = A[(bi + 0) * n + (bk + 0)]; a01 = A[(bi + 0) * n + (bk + 1)]; a02 = A[(bi + 0) * n + (bk + 2)]; a03 = A[(bi + 0) * n + (bk + 3)];
				a10 = A[(bi + 1) * n + (bk + 0)]; a11 = A[(bi + 1) * n + (bk + 1)]; a12 = A[(bi + 1) * n + (bk + 2)]; a13 = A[(bi + 1) * n + (bk + 3)];
				a20 = A[(bi + 2) * n + (bk + 0)]; a21 = A[(bi + 2) * n + (bk + 1)]; a22 = A[(bi + 2) * n + (bk + 2)]; a23 = A[(bi + 2) * n + (bk + 3)];
				a30 = A[(bi + 3) * n + (bk + 0)]; a31 = A[(bi + 3) * n + (bk + 1)]; a32 = A[(bi + 3) * n + (bk + 2)]; a33 = A[(bi + 3) * n + (bk + 3)];

				b00 = B[(bj + 0) * n + (bk + 0)]; b01 = B[(bj + 0) * n + (bk + 1)]; b02 = B[(bj + 0) * n + (bk + 2)]; b03 = B[(bj + 0) * n + (bk + 3)];
				b10 = B[(bj + 1) * n + (bk + 0)]; b11 = B[(bj + 1) * n + (bk + 1)]; b12 = B[(bj + 1) * n + (bk + 2)]; b13 = B[(bj + 1) * n + (bk + 3)];
				b20 = B[(bj + 2) * n + (bk + 0)]; b21 = B[(bj + 2) * n + (bk + 1)]; b22 = B[(bj + 2) * n + (bk + 2)]; b23 = B[(bj + 2) * n + (bk + 3)];
				b30 = B[(bj + 3) * n + (bk + 0)]; b31 = B[(bj + 3) * n + (bk + 1)]; b32 = B[(bj + 3) * n + (bk + 2)]; b33 = B[(bj + 3) * n + (bk + 3)];

				c00 += a00 * b00 + a01 * b01 + a02 * b02 + a03 * b03;
				c01 += a00 * b10 + a01 * b11 + a02 * b12 + a03 * b13;
				c02 += a00 * b20 + a01 * b21 + a02 * b22 + a03 * b23;
				c03 += a00 * b30 + a01 * b31 + a02 * b32 + a03 * b33;

				c10 += a10 * b00 + a11 * b01 + a12 * b02 + a13 * b03;
				c11 += a10 * b10 + a11 * b11 + a12 * b12 + a13 * b13;
				c12 += a10 * b20 + a11 * b21 + a12 * b22 + a13 * b23;
				c13 += a10 * b30 + a11 * b31 + a12 * b32 + a13 * b33;

				c20 += a20 * b00 + a21 * b01 + a22 * b02 + a23 * b03;
				c21 += a20 * b10 + a21 * b11 + a22 * b12 + a23 * b13;
				c22 += a20 * b20 + a21 * b21 + a22 * b22 + a23 * b23;
				c23 += a20 * b30 + a21 * b31 + a22 * b32 + a23 * b33;

				c30 += a30 * b00 + a31 * b01 + a32 * b02 + a33 * b03;
				c31 += a30 * b10 + a31 * b11 + a32 * b12 + a33 * b13;
				c32 += a30 * b20 + a31 * b21 + a32 * b22 + a33 * b23;
				c33 += a30 * b30 + a31 * b31 + a32 * b32 + a33 * b33;
			}

			C[(bi + 0) * n + (bj + 0)] = c00; C[(bi + 0) * n + (bj + 1)] = c01; C[(bi + 0) * n + (bj + 2)] = c02; C[(bi + 0) * n + (bj + 3)] = c03;
			C[(bi + 1) * n + (bj + 0)] = c10; C[(bi + 1) * n + (bj + 1)] = c11; C[(bi + 1) * n + (bj + 2)] = c12; C[(bi + 1) * n + (bj + 3)] = c13;
			C[(bi + 2) * n + (bj + 0)] = c20; C[(bi + 2) * n + (bj + 1)] = c21; C[(bi + 2) * n + (bj + 2)] = c22; C[(bi + 2) * n + (bj + 3)] = c23;
			C[(bi + 3) * n + (bj + 0)] = c30; C[(bi + 3) * n + (bj + 1)] = c31; C[(bi + 3) * n + (bj + 2)] = c32; C[(bi + 3) * n + (bj + 3)] = c33;
		}

	t2 = clock();
	printf("Time for unrolled block multiplication n=%d (block size = %d): %f, C[n/2][n/2] = %f\n", n, bs, (double)(t2 - t1) / CLOCKS_PER_SEC, C[n / 2 * n + n / 2]);
}

void Strassen(double* A, double* B, double* C, int size)
{
	int N = 256; // Размер самого маленького блока 

	if (size <= N)
		blockMultiply(A, B, C, size, false);
	else
	{
		int n = size / 2;

		double* P = new double[7 * n * n];

		double* a = new double[n * n];
		double* b = new double[n * n];

		double* P1, * P2, * P3, * P4, * P5, * P6, * P7;

		P1 = P;
		P2 = P + n * n;
		P3 = P + 2 * n * n;
		P4 = P + 3 * n * n;
		P5 = P + 4 * n * n;
		P6 = P + 5 * n * n;
		P7 = P + 6 * n * n;

		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
			{
				a[i * n + j] = A[i * size + j + n] - A[(i + n) * size + j + n];
				b[i * n + j] = B[(i + n) * size + j] + B[(i + n) * size + j + n];
			}

		Strassen(a, b, P1, n);

		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
			{
				a[i * n + j] = A[i * size + j] + A[(i + n) * size + j + n];
				b[i * n + j] = B[i * size + j] + B[(i + n) * size + j + n];
			}

		Strassen(a, b, P2, n);

		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
			{
				a[i * n + j] = A[i * size + j] - A[(i + n) * size + j];
				b[i * n + j] = B[i * size + j] + B[i * size + j + n];
			}

		Strassen(a, b, P3, n);

		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
			{
				a[i * n + j] = A[i * size + j] + A[i * size + j + n];
				b[i * n + j] = B[(i + n) * size + j + n];
			}

		Strassen(a, b, P4, n);

		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
			{
				a[i * n + j] = A[i * size + j];
				b[i * n + j] = B[i * size + j + n] - B[(i + n) * size + j + n];
			}

		Strassen(a, b, P5, n);

		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
			{
				a[i * n + j] = A[(i + n) * size + j + n];
				b[i * n + j] = B[(i + n) * size + j] - B[i * size + j];
			}

		Strassen(a, b, P6, n);

		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
			{
				a[i * n + j] = A[(i + n) * size + j] + A[(i + n) * size + j + n];
				b[i * n + j] = B[i * size + j];
			}

		Strassen(a, b, P7, n);

		//Собираем C
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
			{
				C[i * size + j] = P1[i * n + j] + P2[i * n + j] - P4[i * n + j] + P6[i * n + j];
				C[i * size + j + n] = P4[i * n + j] + P5[i * n + j];
				C[(i + n) * size + j] = P6[i * n + j] + P7[i * n + j];
				C[(i + n) * size + j + n] = P2[i * n + j] - P3[i * n + j] + P5[i * n + j] - P7[i * n + j];
			}

		delete[] P;
		delete[] a;
		delete[] b;
	}
}

void multiplyIJK(const std::vector<double>& A, const std::vector<double>& B, std::vector<double>& C)
{
	clock_t t1, t2;
	t1 = clock();

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			for (int k = 0; k < n; ++k)
				C[i * n + j] += A[i * n + k] * B[k * n + j];

	t2 = clock();
	printf("Time for n=%d (IJK): %f, C[n/2][n/2] = %f\n", n, (double)(t2 - t1) / CLOCKS_PER_SEC, C[n / 2 * n + n / 2]);
}

void multiplyJKI(const std::vector<double>& A, const std::vector<double>& B, std::vector<double>& C)
{
	clock_t t1, t2;
	t1 = clock();

	for (int j = 0; j < n; ++j)
		for (int k = 0; k < n; ++k)
			for (int i = 0; i < n; ++i)
				C[i * n + j] += A[i * n + k] * B[k * n + j];

	t2 = clock();
	printf("Time for n=%d (JKI): %f, C[n/2][n/2] = %f\n", n, (double)(t2 - t1) / CLOCKS_PER_SEC, C[n / 2 * n + n / 2]);
}

void multiplyKIJ(const std::vector<double>& A, const std::vector<double>& B, std::vector<double>& C)
{
	clock_t t1, t2;
	t1 = clock();

	for (int k = 0; k < n; ++k)
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				C[i * n + j] += A[i * n + k] * B[k * n + j];

	t2 = clock();
	printf("Time for n=%d (KIJ): %f, C[n/2][n/2] = %f\n", n, (double)(t2 - t1) / CLOCKS_PER_SEC, C[n / 2 * n + n / 2]);
}

void multiplyJIK(const std::vector<double>& A, const std::vector<double>& B, std::vector<double>& C)
{
	clock_t t1, t2;
	t1 = clock();

	for (int j = 0; j < n; ++j)
		for (int i = 0; i < n; ++i)		
			for (int k = 0; k < n; ++k)
				C[i * n + j] += A[i * n + k] * B[k * n + j];

	t2 = clock();
	printf("Time for n=%d (JIK): %f, C[n/2][n/2] = %f\n", n, (double)(t2 - t1) / CLOCKS_PER_SEC, C[n / 2 * n + n / 2]);
}

void multiplyIKJ(const std::vector<double>& A, const std::vector<double>& B, std::vector<double>& C)
{
	clock_t t1, t2;
	t1 = clock();

	for (int i = 0; i < n; ++i)
		for (int k = 0; k < n; ++k)
			for (int j = 0; j < n; ++j)			
				C[i * n + j] += A[i * n + k] * B[k * n + j];

	t2 = clock();
	printf("Time for n=%d (IKJ): %f, C[n/2][n/2] = %f\n", n, (double)(t2 - t1) / CLOCKS_PER_SEC, C[n / 2 * n + n / 2]);
}

void multiplyKJI(const std::vector<double>& A, const std::vector<double>& B, std::vector<double>& C)
{
	clock_t t1, t2;
	t1 = clock();

	for (int k = 0; k < n; ++k)
		for (int j = 0; j < n; ++j)
			for (int i = 0; i < n; ++i)			
				C[i * n + j] += A[i * n + k] * B[k * n + j];

	t2 = clock();
	printf("Time for n=%d (KJI): %f, C[n/2][n/2] = %f\n", n, (double)(t2 - t1) / CLOCKS_PER_SEC, C[n / 2 * n + n / 2]);
}

int main()
{
	std::vector<double> A, B, C;
	A.resize(n * n);
	B.resize(n * n);
	C.resize(n * n);

	initMatrix(A, B, C);

/*	multiplyIJK(A, B, C);

	clearC(C);
	multiplyJKI(A, B, C);
	
	clearC(C);
	multiplyKIJ(A, B, C);*/

/*	clearC(C);
	multiplyIKJ(A, B, C);

	clearC(C);
	blockMultiply(A.data(), B.data(), C.data(), n);

	clearC(C);
	blockUnrolledMultiply(A, B, C);

	clearC(C);
	*/
	double t1, t2;
	/*t1 = clock();
	Strassen(A.data(), B.data(), C.data(), n);
	t2 = clock();
	printf("Time for n=%d (Strassen algorithm): %f, C[n/2][n/2] = %f\n", n, (double)(t2 - t1) / CLOCKS_PER_SEC, C[n / 2 * n + n / 2]);
	*/

	t1 = clock();
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1.0,
		A.data(), n, B.data(), n, 0.0, C.data(), n);
	t2 = clock();
	printf("Time for n=%d (MKL algorithm): %f, C[n/2][n/2] = %f\n", n, (double)(t2 - t1) / CLOCKS_PER_SEC, C[n / 2 * n + n / 2]);

	/*clearC(C);
	multiplyKJI(A, B, C);

	clearC(C);
	multiplyJIK(A, B, C);*/

	return EXIT_SUCCESS;
}
