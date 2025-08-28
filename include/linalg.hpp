#pragma once
#include <vector>
#include <stdexcept>
#include <cassert>

namespace linalg {

// ---- Basic utilities (all row-major) ----

// Transpose a row-major matrix S(rows×cols) into T(cols×rows) (row-major).
void transpose_rm(const std::vector<double>& S, int rows, int cols,
                  std::vector<double>& T);

// Matrix multiply: C = A(m×k) * B(k×n)  (row-major everywhere)
void matmul(const std::vector<double>& A, int m, int k,
            const std::vector<double>& B, int k2, int n,
            std::vector<double>& C);

// Matrix-vector: y = A(m×n) * x(n)  (row-major)
void gemv(const std::vector<double>& A, int m, int n,
          const std::vector<double>& x, std::vector<double>& y);

// ---- Reusable LU factorization ----
struct LU {
    int n = 0;                      // matrix size (n×n)
    int info = 0;                   // dgetrf return code
    std::vector<double> data;       // combined L/U factors (row-major)
    std::vector<int>    ipiv;       // pivot indices (1-based in LAPACK docs)

    LU() = default;
    LU(const std::vector<double>& A, int n_) { factor(A, n_); }

    // Factor once (copies A; A is not modified).
    void factor(const std::vector<double>& A, int n_);

    // Solve A x = b (overwrites b with x).
    void solve_inplace(std::vector<double>& b) const;

    // Convenience: return a new vector x, keep input b unchanged.
    std::vector<double> solve(const std::vector<double>& b) const;

    // Multiple right-hand sides.
    // B is n×nrhs, stored row-major with columns being individual RHS vectors:
    //   B(i,j) is at B[i*nrhs + j].
    // Overwrites B with solutions X of the same layout.
    void solve_many_inplace(std::vector<double>& B, int nrhs) const;
};

} // namespace linalg
