#pragma once
#include <cmath>
#include <vector>
#include "solver.hpp"          // father class DGAdvection1D
#include "error_analysis.hpp"  // l2_err2_mass, linf_nodal

// -----------------------------------------------------------------------------
// DGAdvectionDiffusion1D
// PDE: u_t + a u_x = k u_xx on [0,1] with periodic BC
//
// Structure mirrors DGAdvection1D. We reuse:
//   - nodes s, operators Mref, Cref, Mel, LU(Mel)=M_lu, mapping x_of(...)
//   - printing style and error reporting
//
// Diffusion is added via LDG-style auxiliary variable sigma ≈ u_x:
//   1) Mel * sigma = -Cref * u + face terms (alternating flux for q)
//   2) Mel * (du/dt) = Cref * (a u - k sigma) + face terms
// -----------------------------------------------------------------------------
struct DGAdvectionDiffusion1D : public DGAdvection1D {
    double mu;                          // diffusion coefficient
    mutable std::vector<double> sigma; // auxiliary field q ≈ u_x (size: ne*np)

    DGAdvectionDiffusion1D(int p_, int ne_, double a_, double mu_);

    // dudt = RHS(u) (advection + diffusion); same signature as father
    void rhs(const std::vector<double>& u, std::vector<double>& dudt) const;

    // Same signature & defaults policy as father:
    std::vector<double> solve_pde(const double cfl,
                                  const double T,
                                  bool print_results = true,
                                  bool fix_time_step = false,
                                  double dt_fix = 2e-4);

private:
    // Exact solution for sine IC: u(x,0)=sin(2πx) N = 2
    inline double u_exact(double x, double t) const {
        using std::exp;
        using std::sqrt;

        const int N = 2; // periodic images, sufficient here
        const double denom = 1.0 + 400.0 * mu * t;
        const double factor = 1.0 / sqrt(denom);

        double sum = 0.0;
        for (int i = -N; i <= N; ++i) {
            double dx = x - 0.5 - t + i; // shift by periodic image
            sum += exp(-100.0 * dx * dx / denom);
        }

        return factor * sum;
    }
};



