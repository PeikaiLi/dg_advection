#include "quadrature.hpp"
#include <openblas/lapacke.h>
#include <cmath>
#include <vector>
#include <stdexcept>

// Build symmetric tridiagonal J with diag=0, offdiag=b_k = k / sqrt(4k^2-1),
// then eigenvalues -> nodes, weights = 2*(v0)^2 (Golub–Welsch).
void gauss_legendre(int p, std::vector<double>& x, std::vector<double>& w)
{
    const int n = (int)std::ceil((p + 1.0) / 2.0);
    x.assign(n, 0.0);
    w.assign(n, 0.0);

    // d: main diagonal (all zeros), e: off-diagonal (1..n-1)
    std::vector<double> d(n, 0.0), e(n-1, 0.0);
    for (int k = 1; k <= n-1; ++k) e[k-1] = k / std::sqrt(4.0*k*k - 1.0);

    // LAPACKE_dstev: eigenvalues of symmetric tridiagonal + eigenvectors (Z)
    std::vector<double> Z(n*n, 0.0); // column-major
    int info = LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', n, d.data(), e.data(), Z.data(), n);
    // Example: Using LAPACKE_dstev to compute eigenvalues and eigenvectors
    // ---------------------------------------------------------------
    // LAPACKE_dstev solves the symmetric tridiagonal eigenvalue problem.
    // Input:
    //   - d: main diagonal of the tridiagonal matrix (size n)
    //   - e: subdiagonal elements (size n-1)
    //   - Z: workspace for eigenvectors (size n*n, column-major)
    // Parameters:
    //   - jobz = 'V' -> compute both eigenvalues and eigenvectors
    //   - n    = size of the problem
    //   - ldz  = leading dimension of Z (usually n)
    // Output (on return):
    //   - d: overwritten with the eigenvalues (sorted in ascending order)
    //   - Z: contains eigenvectors (each eigenvector is stored as a column)
    //   - info = 0 if success, >0 if algorithm failed to converge
    // ---------------------------------------------------------------

    if (info != 0) throw std::runtime_error("dstev failed: info=" + std::to_string(info));

    // d now holds eigenvalues ascending -> nodes; first row of Z are v0
    for (int i = 0; i < n; ++i) {
        x[i] = d[i];
        double v0 = Z[0 + i*n];        // Z(1,i) in column-major
        w[i] = 2.0 * v0 * v0;
    }
}
