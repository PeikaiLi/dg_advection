#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <numbers>
#include <cassert>
#include "solver.hpp"

#ifndef M_PI
#define M_PI (std::numbers::pi)
#endif

/**
 * Compute Chebyshev–Gauss–Lobatto quadrature weights.
 *
 * Nodes are: x_j = cos(pi * j / p),  j = 0,...,p
 *
 * Quadrature formula:
 *   ∫_{-1}^1 f(x) / sqrt(1-x^2) dx  ≈  Σ w_j f(x_j)
 *
 * In practice (DG on [-1,1]) we use these weights to approximate
 *   ∫_{-1}^1 f(x) dx ≈ Σ w_j f(x_j),
 * because our solution u_h is stored exactly at these nodal points.
 *
 * Formula for weights:
 *   w_0 = w_p = π / (2p)
 *   w_j = π / p   for j=1,...,p-1
 */
std::vector<double> chebyshev_gll_weights(int p) {
    std::vector<double> w(p+1);
    for (int j=0; j<=p; ++j) {
        if (j==0 || j==p) w[j] = M_PI / (2.0 * p);
        else              w[j] = M_PI / p;
    }
    return w;
}

/**
 * Compute discrete L2 error norm squared:
 *
 *   ||err||^2 ≈ Σ_{elements} (h/2) Σ_{j=0}^p w_j * err(x_j^K)^2
 *
 * where:
 *   - err is the error vector, stored nodally (same layout as u)
 *   - h is the element size
 *   - w_j are quadrature weights (depending on nodal choice)
 *
 * Logic: 
 *   In DG with nodal basis, our solution is stored at special nodes
 *   (e.g. Chebyshev–GLL). To approximate integrals we re-use the
 *   same nodes + their quadrature weights. This ensures that:
 *     - The discrete L2 norm corresponds to the continuous L2 inner product
 *     - Error does not "artificially grow" when p increases
 */
double l2_error_squared(const DGAdvection1D& dg,
    const std::vector<double>& err)
{
    std::vector<double> wq = chebyshev_gll_weights(dg.p);
    assert((int)wq.size() == dg.np);
    double acc = 0.0;
    for (int e=0; e<dg.ne; ++e) {
        const double* ve = &err[e*dg.np];  // local error at element e
        for (int j=0; j<dg.np; ++j) {
            acc += (dg.h/2.0) * wq[j] * ve[j] * ve[j];
        }
    }
    return acc;
}


int main() {
    const int    p   = 16;
    const int    ne  = 4;
    const double a   = 1.0;
    const double cfl = 0.35;
    const double T   = 1.0;     // one period for sin(2πx) with a=1

    std::cout << "start..." << std::endl;
    DGAdvection1D dg(p, ne, a);
    
    // initial condition: u0(x) = sin(2πx)
    std::vector<double> u(dg.ne * dg.np, 0.0);
    for (int e = 0; e < dg.ne; ++e)
        for (int j = 0; j < dg.np; ++j) {
            double x = dg.x_of(e, j);
            u[e*dg.np + j] = std::sin(2.0 * std::numbers::pi * x);
        }

    // stable time step for RK4: dt ≈ CFL * h / ((2p+1)|a|)
    double dt = cfl * dg.h / ((2*p + 1) * std::abs(a));
    int nsteps = (int)std::ceil(T / dt);
    dt = T / nsteps; // exactly land on T


    // ERROR 
    for (int s = 0; s < nsteps; ++s) dg.rk4_step(u, dt);

    std::cout << "finished dg..." << std::endl;

    // exact solution: shift by a*T
    std::vector<double> uex(u.size(), 0.0), err(u.size(), 0.0);
    for (int e = 0; e < dg.ne; ++e)
        for (int j = 0; j < dg.np; ++j) {
            double x = dg.x_of(e, j);
            double xe = std::fmod(x - a*T, 1.0);
            if (xe < 0) xe += 1.0;
            uex[e*dg.np + j] = std::sin(2.0 * std::numbers::pi * xe);
            err[e*dg.np + j] = u[e*dg.np + j] - uex[e*dg.np + j];
        }

    double L2 = std::sqrt(l2_error_squared(dg, err));
    std::cout << std::setprecision(6)
              << "p=" << p << ", ne=" << ne
              << ", dt=" << dt
              << ", steps=" << nsteps
              << ",  L2 error: " << L2 << "\n";
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
// L2 error squared: 1.25484e-17
// p=16, ne=4, dt=0.0026455, steps=378,  L2 error: 3.54238e-09