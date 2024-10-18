#include "lu_decomposition.cpp"

int main(){
	int n = 3;
	std::vector<double> A = {2, 3, 1, 4, 7, 1, 6, 18, 5};
	std::vector<double> LU = make_LU_decomposition(A, n);
	test_LU_decomposition(LU, A, n);
}