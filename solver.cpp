#include <iostream>
#include <cassert>
#include <iomanip>
#include <cmath>
#include <openblas/cblas.h>
#include <solver.hpp>
#include "error_analysis.hpp"

DGAdvection1D::DGAdvection1D(int p_, int ne_, double a_)
: p(p_), np(p_+1), ne(ne_), a(a_) {
    h = 1.0 / (ne); 
    // since we have ne = number of elements,
    // assume periodic boundary conditions: ne-th element is same as first element named as e0

    // reference operators
    std::vector<double> gx, gw, phi_Q, dphi_Q;
    build_reference_ops(this->p, this->s, gx, gw, phi_Q, dphi_Q,  this->Mref,  this->Cref);

    // Mel = (h/2) * Mref
    Mel = Mref;
    for (double& v : Mel) v *= 0.5 * h;

    // LU(Mel) for solves in RHS
    M_lu.factor(Mel, np);

    Cel = Cref;
    if (a!=1) {
        for (double& v : Cel) v *= a;
    }
}

void DGAdvection1D::rhs(const std::vector<double>& u, std::vector<double>& dudt) const {
    dudt.assign(np*ne, 0.0);

    // unit basis at faces (nodal basis includes endpoints):
    // left face corresponds to local index 0, right face to np-1
    for (int e = 0; e < ne; ++e) {
        const double* ue = &u[e*np];  // element e solution
        // u is the solution vector for the entire domain
        // if we treat u as a row-major matrix, we have:
        //   u = [u_0, u_1, ..., u_{ne-1}]^T
        //   where u_i is the solution vector for element i
        //   ue is u_e an column vector

        // volume term: r = Cel * ue
        std::vector<double> r(np);  // r is a column vector of length np

        cblas_dgemv(CblasRowMajor, CblasNoTrans,
                    np, np,
                    1.0, Cel.data(), np,
                        ue, 1,
                    0.0, r.data(), 1);


        // upwind numerical flux (periodic BC)
        const int eL = (e - 1 + ne) % ne;
        const int eR = (e + 1) % ne;

        double hat_L, hat_R;
        if (a >= 0.0) {
            hat_L = u[eL*np + (np-1)]; // inflow from left neighbor
            // since it start at 0 np-1 is np-th node.
            hat_R = ue[np-1];         // outflow to right
        } else {
            hat_L = ue[0];           // outflow to left
            hat_R = u[eR*np + 0];     // inflow from right neighbor
        }
        // boundary contribution: -a * (phi(+1)*hat_R - phi(-1)*hat_L)
        // with nodal endpoints: only two entries touched
        r[0]      -= -a * hat_L;
        r[np-1]    -= +a * hat_R;

        // solve Mel * (du/dt)_e = r
        std::vector<double> x = r;   // overwrite with solution
        M_lu.solve_inplace(x);

        // scatter back
        double* de = &dudt[e*np];
        for (int i = 0; i < np; ++i) de[i] = x[i];
    }
}

void DGAdvection1D::rk4_step(std::vector<double>& u, double dt) const {
    std::vector<double> k1(u.size()), k2(u.size()), 
                        k3(u.size()), k4(u.size()), 
                        tmp(u.size());

    rhs(u, k1);

    for (size_t i = 0; i < u.size(); ++i)
        tmp[i] = u[i] + 0.5*dt*k1[i];
    rhs(tmp, k2);

    for (size_t i = 0; i < u.size(); ++i)
        tmp[i] = u[i] + 0.5*dt*k2[i];
    rhs(tmp, k3);

    for (size_t i = 0; i < u.size(); ++i)
        tmp[i] = u[i] + dt*k3[i];
    rhs(tmp, k4);

    for (size_t i = 0; i < u.size(); ++i)
        u[i] += (dt/6.0)*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
}



// return error, numerical solution at T 
std::vector<double> DGAdvection1D::solve_pde(const double cfl, 
    const double T,bool print_results,bool fix_time_step, double dt_fix) {

    std::vector<double> u(ne*np, 0.0);
    err.resize(ne*np);
    uex.resize(ne*np);
    // initialize
    // initial condition: u0(x) = sin(2πx), stored at Chebyshev–GLL nodes
    for (int e = 0; e < ne; ++e)
        for (int j = 0; j < np; ++j) {
            double x = x_of(e, j);
            u[e*np + j] = std::sin(2.0 * std::numbers::pi * x);
        }

    // time-stepping
    // stable time step for RK4: dt ≈ CFL * h / ((2p+1)|a|)
    double dt = cfl * h / ((2*p + 1) * std::abs(a));
    if (fix_time_step) {
        dt = dt_fix;
    }
    int steps = (int)std::ceil(T / dt);
    dt = T / steps; // land exactly at T

    for (int s = 0; s < steps; ++s) rk4_step(u, dt);

    // exact solution
    for (int e = 0; e < ne; ++e)
        for (int j = 0; j < np; ++j) {
            double x = x_of(e, j);
            double xe = std::fmod(x - a*T, 1.0);
            if (xe < 0) xe += 1.0;
            uex[e*np + j] = std::sin(2.0 * std::numbers::pi * xe);
            err[e*np + j] = u[e*np + j] - uex[e*np + j];
        }

    // norms
    double l2_sq = l2_err2_mass(err, ne, np, Mel);
    this->l2    = std::sqrt(l2_sq);
    this->linf  = linf_nodal(err);

    if (print_results) {

        std::cout << std::setprecision(6)
                    << "p=" << p
                    << ", ne=" << ne
                    << ", dt=" << dt
                    << ", steps=" << steps << "\n";

        std::cout << std::scientific << std::setprecision(6);
        std::cout << "  L2 error squared (mass-matrix): " << l2_sq << "\n";
        std::cout << "  L2 error:                       " << this->l2    << "\n";
        std::cout << "  Linf error (nodal):             " << this->linf  << "\n";
    }
    return u;
}