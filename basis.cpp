#include "basis.hpp"
#include <cmath>
#include <numbers>


std::vector<double> chebyshev_nodes(int p) {
    std::vector<double> s(p+1);
    // j = p..0 to match Julia's descending order
    for (int j = p, t = 0; j >= 0; --j, ++t) {
        s[t] = std::cos(std::numbers::pi * j / p);
    }
    return s;
}

void legendre_poly(const std::vector<double>& x,
                   const int p,
                   std::vector<double>& Y,
                   std::vector<double>& dY)
{
    const int nx = (int)x.size();
    Y.assign((p+1)*nx, 0.0);
    dY.assign((p+1)*nx, 0.0);

    // Define a small anonymous function (lambda) called `idx`.
    // It takes (i,k) as input and returns the 1D index in a row-major
    // flattened array of size (nx × (p+1)).
    // The `[&]` means "capture external variables by reference",
    // so the lambda can use `p` from the outer scope.
    // Example: idx(2,1) = 2*(p+1)+1 → position in the vector.
    auto idx = [&](int i, int k){ return i*(p+1) + k; };

    // P0 = 1, P1 = x ; P0' = 0, P1' = 1
    // p0 is the constant polynomial 
    // represented in first column
    // Y[idx(i,0)] = 1.0; that's why we set it to 1
    for (int i = 0; i < nx; ++i) {
        Y[idx(i,0)] = 1.0;
        if (p >= 1) {
            Y[idx(i,1)]  = x[i];
            dY[idx(i,1)] = 1.0;
        }
    }
    // Three-term recurrence:
    // (k+1) P_{k+1} = (2k+1) x*P_k - k*P_{k-1}
    // (k+1) P'_{k+1} = (2k+1)(x*P'_k + P_k) - k*P'_{k-1}
    for (int k = 1; k <= p-1; ++k) {
        const double a = (2.0*k + 1.0) / (k + 1.0);
        const double b = (double)k / (k + 1.0);
        for (int i = 0; i < nx; ++i) {
            Y[idx(i,k+1)]  = a * x[i] * Y[idx(i,k)] - b * Y[idx(i,k-1)];
            dY[idx(i,k+1)] = a * (x[i] * dY[idx(i,k)] + Y[idx(i,k)]) - b * dY[idx(i,k-1)];
        }
    }
}