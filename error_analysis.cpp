#include <numbers>
#include <vector>
#include <cmath>
#include <numeric> // std::accumulate

#include "error_analysis.hpp"
#include "solver.hpp"
#include "solver_advdiff.hpp"

double l2_err2_mass(const std::vector<double>& err, int ne, int np, const std::vector<double>& Mel) {
    double s = 0.0;
    for (int e = 0; e < ne; ++e) {
        const double* v = &err[e*np];
        for (int i = 0; i < np; ++i)
            for (int j = 0; j < np; ++j)
                s += v[i] * Mel[i*np + j] * v[j];
    }
    return s;
}

double linf_nodal(const std::vector<double>& err) {
    double m = 0.0; for (double v : err) m = std::max(m, std::abs(v)); return m;
}


Metrics run_once_DGAdvection1D(int p, int ne, double a, double T, double cfl, bool fix_time_step, double dt_fix) {
    DGAdvection1D dg(p, ne, a);
    dg.solve_pde(cfl, T,false, fix_time_step, dt_fix);
    Metrics m;
    m.L2   = dg.l2;
    m.Linf = dg.linf;
    return m;
}

Metrics run_once_DGAdvection1D_Diff(int p, int ne, double mu,double a, double T, double cfl, bool fix_time_step, double dt_fix) {
    DGAdvectionDiffusion1D dg(p, ne, a, mu);
    dg.solve_pde(cfl, T,false, fix_time_step, dt_fix);
    Metrics m;
    m.L2   = dg.l2;
    m.Linf = dg.linf;
    return m;
}


// ---------- least-squares slope y = m x + b ----------
double slope_ls(const std::vector<double>& x, const std::vector<double>& y) {
    const int n = (int)x.size();
    double mx = std::accumulate(x.begin(), x.end(), 0.0) / n;
    double my = std::accumulate(y.begin(), y.end(), 0.0) / n;
    double num = 0.0, den = 0.0;
    for (int i = 0; i < n; ++i) { double dx = x[i]-mx, dy = y[i]-my; num += dx*dy; den += dx*dx; }
    return num/den;
}

// slope between the two neighboring points with maximum slope (by absolute value)
double slope_max_adjacent(const std::vector<double>& x,
                                 const std::vector<double>& y) {
    if (x.size() < 2) return 0.0;
    double maxSlope = (y[1]-y[0]) / (x[1]-x[0]);
    for (size_t i = 1; i+1 < x.size(); ++i) {
        double s = (y[i+1]-y[i]) / (x[i+1]-x[i]);
        if (std::abs(s) > std::abs(maxSlope))
            maxSlope = s;
    }
    return maxSlope;
}
