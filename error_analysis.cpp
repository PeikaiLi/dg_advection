#include <numbers>
#include <vector>
#include <cmath>
#include "error_analysis.hpp"
#include "solver.hpp"

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


Metrics run_once(int p, int ne, double a, double T, double cfl, bool fix_time_step, double dt_fix) {
    DGAdvection1D dg(p, ne, a);
    dg.solve_pde(cfl, T,false, fix_time_step, dt_fix);
    Metrics m;
    m.L2   = dg.l2;
    m.Linf = dg.linf;
    return m;
}
