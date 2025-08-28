#include <iostream>
#include <vector>
#include <cmath>
#include "basis.hpp"
#include "quadrature.hpp"

int main() {
    int p = 5;
    // 1) 测 Legendre 递推是否稳定
    std::vector<double> s = chebyshev_nodes(p);
    std::vector<double> Y, dY;
    legendre_poly(s, p, Y, dY);

    auto idx = [&](int i, int k){ return i*(p+1)+k; };// row first
    std::cout << "s[0]=" << s[0] << " (expect -1)\n";
    std::cout << "P0(s[0])=" << Y[idx(0,0)] << " (expect ~1)\n";
    std::cout << "P1(s[0])=" << Y[idx(0,1)] << " (expect ~s[0])\n";

    // 2) Gauss–Legendre check ∫ P_i*P_j == 2/(2i+1) δ_ij
    std::vector<double> gx, gw;
    gauss_legendre(2*p, gx, gw); // degree 2p -> exactly integrates P_i P_j up to i+j<=2p

    double max_err = 0.0;
    for (int i = 0; i <= p; ++i) {
        for (int j = 0; j <= p; ++j) {
            // evaluate P_i, P_j at gx
            std::vector<double> Yg, dYg;
            legendre_poly(gx, std::max(i,j), Yg, dYg);
            double I = 0.0;
            for (int m = 0; m < (int)gx.size(); ++m) {
                double Pi = Yg[m*(std::max(i,j)+1)+i];
                double Pj = Yg[m*(std::max(i,j)+1)+j];
                I += gw[m] * Pi * Pj;
            }
            double exact = (i==j) ? 2.0/(2*i+1) : 0.0;
            max_err = std::max(max_err, std::abs(I - exact));
        }
    }
    std::cout << "Max orthogonality error = " << max_err << "\n";
    return 0;
}
// g++ -std=c++20 example/check_basis.cpp basis.cpp quadrature.cpp -Iinclude -lopenblas -o check_basis.exe

// ./check_basis.exe