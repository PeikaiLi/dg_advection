#include "solver_advdiff.hpp"
#include <cmath>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <openblas/cblas.h>

// ctor: reuse father's setup; just allocate sigma
DGAdvectionDiffusion1D::DGAdvectionDiffusion1D(int p_, int ne_, double a_, double mu_)
: DGAdvection1D(p_, ne_, a_), mu(mu_) {
    sigma.assign(ne * np, 0.0);
}

// -----------------------------------------------------------------------------
// RHS(u):
// Step 1: solve for sigma (q ≈ u_x)
//   rhs_q = -Cel * u  + faces (alternating flux):
//     left-face  contrib:  -u_e(0)
//     right-face contrib:  +u_{e+1}(0)
//   sigma = Mel^{-1} * rhs_q   (use father's LU)
//
// Step 2: residual r for u-equation:
//   r = Cel * (a*u - mu*sigma)  + faces
//   faces:
//     - advection on u: upwind (same as father)
//     - diffusion on q: alternating  (left uses neighbor-right, right uses this-right)
//
// Strong-form nodal injection (GLL endpoints): only r[0], r[np-1] touched.
// -----------------------------------------------------------------------------
void DGAdvectionDiffusion1D::rhs(const std::vector<double>& u, std::vector<double>& dudt) const {
    dudt.assign(ne * np, 0.0);

    // ---- Step 1: sigma ----
    for (int e = 0; e < ne; ++e) {
        const double* ue = &u[e * np];
        std::vector<double> rhs_q(np, 0.0);

        // rhs_q = -Cel * ue
        cblas_dgemv(CblasRowMajor, CblasNoTrans,
                    np, np,
                    -1.0, Cel.data(), np,
                    ue, 1,
                    0.0, rhs_q.data(), 1);

        const int eR = (e + 1) % ne;

        // alternating flux for q:
        rhs_q[0]      += -ue[0];          // left face: -u_e(left)
        rhs_q[np - 1] += u[eR * np + 0]; // right face: +u_{eR}(left)

        // sigma_e = Mel^{-1} rhs_q
        std::vector<double> qloc = rhs_q;
        M_lu.solve_inplace(qloc);
        for (int i = 0; i < np; ++i) sigma[e * np + i] = qloc[i];
    }

    // ---- Step 2: r for u ----
    for (int e = 0; e < ne; ++e) {
        const int eL = (e - 1 + ne) % ne;
        const int eR = (e + 1) % ne;

        const double* ue = &u[e * np];
        const double* sigma_e = &sigma[e * np];

        // f = a * u - mu * sigma
        std::vector<double> f(np);// a = 1 here
        for (int i = 0; i < np; ++i) f[i] = a * ue[i] - mu * sigma_e[i];

        // r = Cel * f
        std::vector<double> r(np, 0.0);
        cblas_dgemv(CblasRowMajor, CblasNoTrans,
                    np, np,
                    1.0, Cel.data(), np,
                    f.data(), 1,
                    0.0, r.data(), 1);

        // face fluxes (strong injection on endpoints)
        // here a = 1, just want to little bit general
        // upwind for advection(u):
        const double u_L = (a >= 0.0) ? u[eL * np + (np - 1)] : ue[0];
        const double u_R = (a >= 0.0) ? ue[np - 1]            : u[eR * np + 0];

        // alternating for diffusion(q):
        const double sigma_e_L = sigma[eL * np + (np - 1)]; // neighbor's right
        const double sigma_e_R = sigma[e   * np + (np - 1)]; // this element's right

        // injection signs: r[0] += a*u*_L - mu*q*_L;  r[np-1] -= a*u*_R - mu*q*_R
        r[0]      += a * u_L - mu * sigma_e_L;
        r[np - 1] -= a * u_R - mu * sigma_e_R;

        // (du/dt)_e = Mel^{-1} r_e
        std::vector<double> du_e = r;
        M_lu.solve_inplace(du_e);

        double* de = &dudt[e * np];
        for (int i = 0; i < np; ++i) de[i] = du_e[i];
    }
}


// Same user-facing behavior as father; only dt logic adds diffusion stability
std::vector<double> DGAdvectionDiffusion1D::solve_pde(const double cfl,
                                                      const double T,
                                                      bool print_results,
                                                      bool fix_time_step,
                                                      double dt_fix) {
    std::vector<double> u(ne * np, 0.0);
    err.resize(ne * np);
    uex.resize(ne * np);

    // initial condition: u0(x) = exp(-(x - 0.5)^2 / 0.1^2)
    for (int e = 0; e < ne; ++e)
        for (int j = 0; j < np; ++j) {
            const double x = x_of(e, j);
            u[e * np + j] = std::exp(-(x - 0.5)*(x - 0.5) / (0.1*0.1));
        }

    // dt from advection & diffusion (explicit)
    double dt;
    if (fix_time_step) {
        dt = dt_fix;
    } else {
        const double eps = 1e-300;
        const double adv = (std::abs(a) > eps) ? (h / ((2 * p + 1) * std::abs(a))) : 1e300;
        const double dif = (mu > eps) ? (0.5 * h * h / (mu * (2 * p + 1) * (2 * p + 1))) : 1e300;
        dt = cfl * std::min(adv, dif);
    }
    int steps = (int)std::ceil(T / dt);
    dt = T / steps;

    for (int s = 0; s < steps; ++s) rk4_step(u, dt);

    // exact solution and error
    for (int e = 0; e < ne; ++e)
        for (int j = 0; j < np; ++j) {
            const double x = x_of(e, j);
            uex[e * np + j] = u_exact(x, T);
            err[e * np + j] = u[e * np + j] - uex[e * np + j];
        }

    // norms (same as father)
    const double l2_sq = l2_err2_mass(err, ne, np, Mel);
    this->l2   = std::sqrt(l2_sq);
    this->linf = linf_nodal(err);

    if (print_results) {
        std::cout << std::setprecision(6)
                  << "p=" << p
                  << ", ne=" << ne
                  << ", dt=" << (T / steps)
                  << ", steps=" << steps << "\n";

        std::cout << std::scientific << std::setprecision(6);
        std::cout << "  L2 error squared (mass-matrix): " << l2_sq   << "\n";
        std::cout << "  L2 error:                       " << this->l2   << "\n";
        std::cout << "  Linf error (nodal):             " << this->linf << "\n";
    }
    return u;
}
