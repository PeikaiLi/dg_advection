#include "element.hpp"
#include "basis.hpp"
#include "quadrature.hpp"
#include "linalg.hpp"            // only transpose_rm
#include <openblas/lapacke.h>
#include <openblas/cblas.h>
#include <cassert>
#include <vector>
#include <algorithm>

using std::vector;

void build_reference_ops(int p,
                         vector<double>& s,
                         vector<double>& gx, vector<double>& gw,
                         vector<double>& phi_Q, vector<double>& dphi_Q,
                         vector<double>& Mref,  vector<double>& Cref)
{
    const int n = p + 1;

    // 1) Nodes
    s = chebyshev_nodes(p);          // Chebyshev nodes, size = n (descending order)
    gauss_legendre(2*p, gx, gw);     // Gauss–Legendre quadrature
    const int ng = (int)gx.size();
    assert((int)gw.size() == ng);

    // 2) Legendre basis table
    vector<double> Ys, dYs, Yg, dYg;
    legendre_poly(s,  p, Ys, dYs);   // Ys, dYs: (n x n), row-major
    legendre_poly(gx, p, Yg, dYg);   // Yg, dYg: (ng x n), row-major
    // Ys,Yg is value of Legendre polynomials at Chebyshev and Gauss-Legendre nodes
    // dYs,dYg is value of derivatives of Legendre polynomials at Chebyshev and Gauss-Legendre nodes

    // 3) Right division: phi_Q = Yg / Ys ; dphi_Q = dYg / Ys
    // phi_Q = Yg / Ys equivalent to: Yg = phi_Q * Ys
    // equivalent to Ys^T *  phi_Q^T =  Yg^T
    // Strategy: LU factorization of (Ys^T), then two solves with dgetrs(Solve BT * X = AT)
    vector<double> BT;                         // BT = Ys^T (row-major)
    linalg::transpose_rm(Ys, n, n, BT);

    std::vector<lapack_int> ipiv(n);
    lapack_int N = (lapack_int)n;
    lapack_int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, N, N, BT.data(), N, ipiv.data());
    assert(info == 0 && "LU(Ys^T) failed");

    // First solve: phi_Q = Yg / Ys
    vector<double> AT;                          // AT = (Yg)^T : (n x ng)
    linalg::transpose_rm(Yg, ng, n, AT);
    // LAPACK row-major convention: ldb = nrhs (here nrhs = ng)
    info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', N, (lapack_int)ng,
                          BT.data(), N, ipiv.data(),
                          AT.data(), (lapack_int)ng);// AT will be overwritten with phi_Q
    assert(info == 0 && "solve for phi failed");
    linalg::transpose_rm(AT, n, ng, phi_Q);    // back to row-major  is (ng x n)

    // Second solve: dphi_Q = dYg / Ys
    vector<double> AT2;
    linalg::transpose_rm(dYg, ng, n, AT2);     // (n x ng)
    info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', N, (lapack_int)ng,
                          BT.data(), N, ipiv.data(),
                          AT2.data(), (lapack_int)ng);
    assert(info == 0 && "solve for dphi failed");
    linalg::transpose_rm(AT2, n, ng, dphi_Q);  // dphi_Q is (ng x n)

    // 4) Reference matrices (using BLAS for performance)
    //    Mref = Phi^T * diag(w) * Phi
    //    Cref = dPhi^T * diag(w) * Phi
    Mref.assign(n*n, 0.0);
    Cref.assign(n*n, 0.0);

    // Scale rows by quadrature weights: PhiW = diag(w) * Phi ; same for dPhiW
    vector<double> PhiW  = phi_Q;   // (ng x n)
    vector<double> dPhiW = dphi_Q;  // (ng x n)
    for (int m = 0; m < ng; ++m) {
        const double w = gw[m];
        double* r1 = &PhiW [m*n];
        double* r2 = &dPhiW[m*n];
        for (int j = 0; j < n; ++j) { r1[j] *= w; r2[j] *= w; }
    }

    // Compute Mref = phi_Q^T * PhiW using BLAS dgemm
    cblas_dgemm(
        CblasRowMajor,    // Layout: matrices are stored row-major (C style)
        CblasTrans,       // TransA: use transpose of A (phi_Q^T)
        CblasNoTrans,     // TransB: use B (PhiW) as is
        n,                // M: number of rows of result C (n x n)
        n,                // N: number of cols of result C (n x n)
        ng,               // K: inner dimension = cols of op(A) = rows of op(B)
        1.0,              // alpha: scalar multiplier for A*B
        phi_Q.data(), n,  // A: pointer to phi_Q (ng x n), lda = row length = n
        PhiW.data(), n,   // B: pointer to PhiW (ng x n), ldb = row length = n
        0.0,              // beta: multiply existing C by beta (0.0 = overwrite)
        Mref.data(), n    // C: result matrix (n x n), ldc = row length = n
    );


    // Use GEMM: Cref = dphi_Q^T * PhiW
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                n, n, ng, 1.0,
                dphi_Q.data(), n,
                PhiW.data(),   n,
                0.0,
                Cref.data(),   n);

    // Note (assembly on physical element):
    //   Mel = (h/2) * Mref
    //   Cel = Cref
}


// g++ -std=c++20 \
//     basis.cpp quadrature.cpp element.cpp linalg.cpp example/check_basis.cpp \
//     -Iinclude -lopenblas -o check_basis.exe
