#include "linalg.hpp"
#include <openblas/cblas.h>
#include <openblas/lapacke.h>

namespace linalg {

void transpose_rm(const std::vector<double>& S, int rows, int cols,
                  std::vector<double>& T) {
    assert((int)S.size() == rows * cols);
    T.assign(cols * rows, 0.0);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            T[j*rows + i] = S[i*cols + j];
}

void matmul(const std::vector<double>& A, int m, int k,
            const std::vector<double>& B, int k2, int n,
            std::vector<double>& C) {
    assert(k == k2);
    assert((int)A.size() == m*k);
    assert((int)B.size() == k*n);
    C.assign(m * n, 0.0);

    // Row-major leading dimensions:
    //   A(m×k): lda = k
    //   B(k×n): ldb = n
    //   C(m×n): ldc = n
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, k,
                1.0, A.data(), k,
                     B.data(), n,
                0.0, C.data(), n);
}

void gemv(const std::vector<double>& A, int m, int n,
          const std::vector<double>& x, std::vector<double>& y) {
    assert((int)A.size() == m*n);
    assert((int)x.size() == n);
    y.assign(m, 0.0);

    // Row-major leading dimension for A(m×n) is n.
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
                m, n,
                1.0, A.data(), n,
                x.data(), 1,
                0.0, y.data(), 1);      
}

// ---------------- LU ----------------

void LU::factor(const std::vector<double>& A, int n_) {
    n = n_;
    assert((int)A.size() == n*n);
    data = A;                 // deep copy; keep A intact
    ipiv.assign(n, 0);

    // Row-major LU factorization with partial pivoting.
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n,
                          data.data(), n, ipiv.data());

    if (info < 0) {
        throw std::invalid_argument(
            "dgetrf: argument " + std::to_string(-info) + " had an illegal value");
    }

    if (info > 0) {
        throw std::runtime_error(
            "LU is singular: zero pivot at index " + std::to_string(info));
    }
}

void LU::solve_inplace(std::vector<double>& b) const {

    assert((int)b.size() == n);
    const int nrhs = 1;  
    // number of RHS columns; here we solve for a single vector

    // Solve A * X = B using the LU factors produced by dgetrf.
    //
    // LAPACKE_dgetrs(
    //   order,   // Memory layout: LAPACK_ROW_MAJOR because our arrays are row-major
    //   trans,   // 'N' -> solve A * X = B
    //            // 'T' -> solve A^T * X = B
    //            // 'C' -> solve A^H * X = B (complex only; same as 'T' for real)
    //   n,       // Order of A (A is n-by-n)
    //   nrhs,    // Number of right-hand sides (columns of B / X)
    //   a,       // Pointer to LU factors from dgetrf (combined L and U, unit diag on L)
    //   lda,     // Leading dimension of A(A is n-by-n so it is n)
    //   ipiv,    // Pivot indices from dgetrf (length n, 1-based)
    //   b,       // On entry: B (n-by-nrhs). On exit: X (solution), overwriting B.
    //   ldb      // Because we are using row-major,
                  // the leading dimension of B is the row stride = number of columns of B, i.e. nrhs (not n).
    // )
    //
    // IMPORTANT row-major gotcha:
    //  - For row-major layout, the "leading dimension" is the LENGTH OF THE LAST INDEX.
    //    So for B (size n-by-nrhs), ldb must be nrhs (NOT n).
    //  - For column-major (Fortran) layout, it would be ldb = n.

    // Row-major rule for GETRS:
    //   ldb = nrhs  (NOT n)
    int info_s = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N',
                                n, nrhs,
                                data.data(), n, ipiv.data(),
                                b.data(), nrhs);
    if (info_s != 0)
        throw std::runtime_error("dgetrs failed: info=" + std::to_string(info_s));
}

std::vector<double> LU::solve(const std::vector<double>& b) const {
    std::vector<double> x = b;
    solve_inplace(x);
    return x;
}

void LU::solve_many_inplace(std::vector<double>& B, int nrhs) const {
    assert((int)B.size() == n * nrhs);
    // In row-major storage, B is n×nrhs with ldb = nrhs,
    // and each column j is the j-th RHS stacked every 'nrhs' entries.
    int info_s = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N',
                                n, nrhs,
                                data.data(), n, ipiv.data(),
                                B.data(), nrhs);
    if (info_s != 0)
        throw std::runtime_error("dgetrs(many) failed: info=" + std::to_string(info_s));
}

} // namespace linalg
