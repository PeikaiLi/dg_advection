#include <iostream>
#include <openblas/cblas.h>
#include <openblas/lapacke.h>

int main() {
    int n = 2;
    double A[4] = {4, 1, 2, 3}; // column-major 2x2
    int ipiv[2];
    int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, A, n, ipiv);

    if (info == 0) {
        std::cout << "LU factorization success!\n";
        std::cout << "Matrix A (LU combined):\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                std::cout << A[i + j*n] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "Pivot indices: " << ipiv[0] << " " << ipiv[1] << "\n";
    } else {
        std::cout << "LU failed, info=" << info << "\n";
    }
    return 0;
}

// g++ -O2 -std=c++17 lapack_test.cpp -lopenblas -o test.exe
// make // complie
// make run
// make clean