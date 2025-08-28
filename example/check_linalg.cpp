#include "linalg.hpp"
#include <iostream>
#include <iomanip>

using namespace linalg;

int main() {
    // A (3×3), row-major
    std::vector<double> A = {
         2,  1, 1,
         4, -6, 0,
        -2,  7, 2
    };
    const int n = 3;

    // 1) Factor once
    LU A_lu(A, n);
    // 2a) Single RHS: solve A x = b1
    std::vector<double> b1 = {5, -2, 9};
    auto x1 = A_lu.solve(b1);
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "x1 = [" << x1[0] << ", " << x1[1] << ", " << x1[2] << "]\n";

    // 2b) Multiple RHS at once: columns are b1 and b2
    const int nrhs = 2;
    std::vector<double> B(n * nrhs);
    std::vector<double> b2 = {1, 2, 3};
    for (int i = 0; i < n; ++i) {
        B[i*nrhs + 0] = b1[i]; // column 0 = b1
        B[i*nrhs + 1] = b2[i]; // column 1 = b2
    }
    A_lu.solve_many_inplace(B, nrhs);
    std::cout << "X(col 0) = [" << B[0*nrhs+0] << ", " << B[1*nrhs+0] << ", " << B[2*nrhs+0] << "]\n";
    std::cout << "X(col 1) = [" << B[0*nrhs+1] << ", " << B[1*nrhs+1] << ", " << B[2*nrhs+1] << "]\n";

    // 3) Quick check: A * x1 ≈ b1
    std::vector<double> Ax;
    matmul(A, n, n, x1, n, 1, Ax); // treat x1 as n×1, returns Ax as n×1
    std::cout << "A*x1 = [" << Ax[0] << ", " << Ax[1] << ", " << Ax[2] << "]  (should match b1)\n";

    return 0;
}
//  g++ -std=c++20 example/check_linalg.cpp linalg.cpp -Iinclude -lopenblas -o check_linalg.exe
// ./check_linalg.exe
//  basis.cpp quadrature.cpp // do not need

