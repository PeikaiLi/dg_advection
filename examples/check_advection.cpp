#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <numbers>
#include "solver.hpp"

/*
  L2 error via the mass matrix (what/why):

  In each element K we represent the DG function in a nodal Lagrange basis
      u_h(r) = sum_{j=0}^p u_j φ_j(r),  r ∈ [-1,1].
  The continuous L2 inner-product in this finite-dimensional space is represented
  EXACTLY by the (element) mass matrix:
      (M_ref)_{ij} = ∫_{-1}^1 φ_i φ_j dr,   M_el = (h/2) * M_ref.
  Therefore, for two coefficient vectors v,w (i.e., nodal values),
      ∫_K v_h w_h dx  =  v^T M_el w.
  Taking v = w = e (the nodal error/interpolant), we get
      ||e_h||_{L2(K)}^2 = e^T M_el e.
  This matches the DG inner product by construction and is mathematically exact
  for polynomials up to degree 2p (since M_ref was built with Gauss–Legendre).
*/

// L2(error)^2 = sum_e e_e^T * Mel * e_e
static double l2_error_squared_mass(const DGAdvection1D& dg,
                                    const std::vector<double>& err) {
    double acc = 0.0;
    for (int e = 0; e < dg.ne; ++e) {
        const double* ve = &err[e*dg.np];
        for (int i = 0; i < dg.np; ++i)
            for (int j = 0; j < dg.np; ++j)
                acc += ve[i] * dg.Mel[i*dg.np + j] * ve[j];
    }
    return acc;
}

// Nodal max-norm: max over stored DOF locations (Chebyshev–GLL points).
// This is not the continuous sup-norm on [0,1], but it’s the standard DG
// "nodal L∞" and is what you can compute directly from the DOFs.
static double linf_error_nodal(const DGAdvection1D& dg,
                               const std::vector<double>& err) {
    double m = 0.0;
    for (double v : err) m = std::max(m, std::abs(v));
    return m;
}

int main() {
    const int    p   = 16;
    const int    Ne  = 16;     // number of elements
    const double a   = 1.0;   // advection speed
    const double cfl = 0.35;
    const double T   = 1.0;   // one period for sin(2πx) with a=1

    std::cout << "start...\n";
    DGAdvection1D dg(p, Ne, a);

    // initial condition: u0(x) = sin(2πx), stored at Chebyshev–GLL nodes
    std::vector<double> u(dg.ne * dg.np, 0.0);
    for (int e = 0; e < dg.ne; ++e)
        for (int j = 0; j < dg.np; ++j) {
            double x = dg.x_of(e, j);
            u[e*dg.np + j] = std::sin(2.0 * std::numbers::pi * x);
        }

    // stable time step for RK4: dt ≈ CFL * h / ((2p+1)|a|)
    double dt    = cfl * dg.h / ((2*p + 1) * std::abs(a));
    int    steps = (int)std::ceil(T / dt);
    dt = T / steps; // land exactly at T

    for (int s = 0; s < steps; ++s) dg.rk4_step(u, dt);
    std::cout << "finished dg...\n";

    // exact solution at final time: shift by a*T (periodic domain [0,1])
    std::vector<double> uex(u.size()), err(u.size());
    for (int e = 0; e < dg.ne; ++e)
        for (int j = 0; j < dg.np; ++j) {
            double x  = dg.x_of(e, j);
            double xe = std::fmod(x - a*T, 1.0);  // wrap to [0,1)
            if (xe < 0) xe += 1.0;
            uex[e*dg.np + j] = std::sin(2.0 * std::numbers::pi * xe);
            err[e*dg.np + j] = u[e*dg.np + j] - uex[e*dg.np + j];
        }

    // norms
    double L2_sq = l2_error_squared_mass(dg, err);
    double L2    = std::sqrt(L2_sq);
    double Linf  = linf_error_nodal(dg, err);

    std::cout << std::setprecision(6)
              << "p=" << p
              << ", ne=" << Ne
              << ", dt=" << dt
              << ", steps=" << steps << "\n";
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "  L2 error squared (mass-matrix): " << L2_sq << "\n";
    std::cout << "  L2 error:                       " << L2    << "\n";
    std::cout << "  Linf error (nodal):             " << Linf  << "\n";
    return 0;
}

// g++ -std=c++20 \
//   basis.cpp quadrature.cpp element.cpp linalg.cpp solver.cpp \
//   examples/check_advection.cpp \
//   -Iinclude -lopenblas -o check_advection.exe
// ./check_advection.exe
// $ ./check_advection.exe
// start...
// finished dg...
// p=16, ne=4, dt=0.0026455, steps=378
//   L2 error squared (mass-matrix): 7.988583e-18
//   L2 error:                       2.826408e-09
//   Linf error (nodal):             3.997140e-09

// p=16, ne=4, dt=0.0026455, steps=378
//   L2 error squared (mass-matrix): 7.988583e-18
//   L2 error:                       2.826408e-09
//   Linf error (nodal):             3.997140e-09