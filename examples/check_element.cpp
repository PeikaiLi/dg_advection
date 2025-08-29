#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include "element.hpp"   // build_reference_ops

// max |A|
static double max_abs(const std::vector<double>& A) {
    double m = 0.0;
    for (double v : A) m = std::max(m, std::abs(v));
    return m;
}

// max |A - B|, A,B are row-major n×n
static double max_abs_diff(const std::vector<double>& A,
                           const std::vector<double>& B, int n) {
    double m = 0.0;
    for (int i = 0; i < n*n; ++i) m = std::max(m, std::abs(A[i] - B[i]));
    return m;
}

// C = A^T (row-major)
static void transpose_rm(const std::vector<double>& A, int rows, int cols,
                         std::vector<double>& C) {
    C.assign(cols*rows, 0.0);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            C[j*rows + i] = A[i*cols + j];
}

int main() {
    std::cout << std::setprecision(12);
    for (int p : {2, 3, 5}) {
        const int n = p + 1;
        std::vector<double> s, gx, gw, phi_Q, dphi_Q, Mref, Cref;
        build_reference_ops(p, s, gx, gw, phi_Q, dphi_Q, Mref, Cref);

        // 1) symmetric ‖M - M^T‖∞
        std::vector<double> MT;
        transpose_rm(Mref, n, n, MT);
        double sym_err = max_abs_diff(Mref, MT, n);

        // 2) integral by parts C + C^T = S
        //    Chebyshev-G-Lobatto only has 2 non-zero entries:
        //    located on the end points.
        //    S(n-1,n-1) = +1,  S(0,0) = -1
        std::vector<double> CT;
        transpose_rm(Cref, n, n, CT);

        std::vector<double> S(n*n, 0.0);
        S[(n-1)*n + (n-1)] = +1.0;  // +1 at x=+1
        S[0]               = -1.0;  // -1 at x=-1

        // E = (C + C^T) - S
        std::vector<double> E(n*n, 0.0);
        for (int i = 0; i < n*n; ++i) E[i] = Cref[i] + CT[i] - S[i];
        double ibp_err = max_abs(E);   // integration-by-parts error


        std::cout << "p = " << p
                  << "  n = " << n
                  << "  ng = " << (int)gx.size() << "\n";
        std::cout << "  ||M - M^T||_inf = " << sym_err << "\n";
        std::cout << "  ||(C + C^T) - S||_inf = " << ibp_err << "\n\n";
    }

    return 0;
}
// g++ -std=c++20 \
//   basis.cpp quadrature.cpp element.cpp linalg.cpp examples/check_element.cpp \
//   -Iinclude -lopenblas -o check_element.exe
// ./check_element.exe


// p = 2  n = 3  ng = 3
//   ||M - M^T||_inf = 2.77555756156e-17
//   ||(C + C^T) - S||_inf = 6.66133814775e-16

// p = 3  n = 4  ng = 4
//   ||M - M^T||_inf = 2.25514051877e-17
//   ||(C + C^T) - S||_inf = 1.00740083074e-15

// p = 5  n = 6  ng = 6
//   ||M - M^T||_inf = 1.38777878078e-17
//   ||(C + C^T) - S||_inf = 4.4408920985e-15