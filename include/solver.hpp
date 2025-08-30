#pragma once
#include <vector>
#include "element.hpp"
#include "linalg.hpp"
#include "basis.hpp"

struct DGAdvection1D {
    int p, np, ne;         // degree, dofs/elem(degrees of freedom per element), number of elements
    double a, h;           // wave speed, element size
    // reference operators
    std::vector<double> Mref, Cref;
    // element mass & its LU for solves
    std::vector<double> Mel;// M on element
    linalg::LU M_lu;
    // volume operator Cel = a * Cref //  C on element
    std::vector<double> Cel;
    // nodal points for mapping to x (desc: -1..1, Chebyshev-G-Lobatto)
    std::vector<double> s;

    DGAdvection1D(int p_, int np_, double a_);

    // dudt = RHS(u)
    virtual void rhs(const std::vector<double>& u, std::vector<double>& dudt) const;

    // single RK4 step: u <- u + dt * Φ(u)
    void rk4_step(std::vector<double>& u, double dt) const;

    // physical coordinate of element e, local node j
    inline double x_of(int e, int j) const {
        // r in [-1,1] -> x in [e*h,(e+1)h]
        return e*h + 0.5*h*(s[j] + 1.0);
    }

    std::vector<double> err;
    std::vector<double> uex;
    double l2, linf;
    // return error, numerical solution at T 
    std::vector<double> solve_pde(const double cfl, 
        const double T,bool print_results = true, bool fix_time_step = false, double dt_fix = 2*1e-4);
};
