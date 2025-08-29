#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <numbers>
#include <numeric>
#include "solver.hpp"

// ---------- norms (L2 via mass matrix; nodal Linf) ----------
static double l2_err2_mass(const DGAdvection1D& dg, const std::vector<double>& err) {
    double s = 0.0;
    for (int e = 0; e < dg.ne; ++e) {
        const double* v = &err[e*dg.np];
        for (int i = 0; i < dg.np; ++i)
            for (int j = 0; j < dg.np; ++j)
                s += v[i] * dg.Mel[i*dg.np + j] * v[j];
    }
    return s;
}
static double linf_nodal(const std::vector<double>& err) {
    double m = 0.0; for (double v : err) m = std::max(m, std::abs(v)); return m;
}

// ---------- one run ----------
struct Metrics { double L2; double Linf; };
static Metrics run_once(int p, int ne, double a, double T, double cfl) {
    DGAdvection1D dg(p, ne, a);

    std::vector<double> u(dg.ne*dg.np);
    for (int e = 0; e < dg.ne; ++e)
        for (int j = 0; j < dg.np; ++j) {
            double x = dg.x_of(e, j);
            u[e*dg.np + j] = std::sin(2.0 * std::numbers::pi * x);
        }

    double dt = cfl * dg.h / ((2*p + 1) * std::abs(a));
    int steps = (int)std::ceil(T / dt);
    dt = T / steps;
    for (int s = 0; s < steps; ++s) dg.rk4_step(u, dt);

    std::vector<double> err(u.size());
    for (int e = 0; e < dg.ne; ++e)
        for (int j = 0; j < dg.np; ++j) {
            double x  = dg.x_of(e, j);
            double xe = std::fmod(x - a*T, 1.0); if (xe < 0) xe += 1.0;
            double ue = std::sin(2.0 * std::numbers::pi * xe);
            err[e*dg.np + j] = u[e*dg.np + j] - ue;
        }

    Metrics m;
    m.L2   = std::sqrt(l2_err2_mass(dg, err));
    m.Linf = linf_nodal(err);
    return m;
}

// ---------- least-squares slope y = m x + b ----------
static double slope_ls(const std::vector<double>& x, const std::vector<double>& y) {
    const int n = (int)x.size();
    double mx = std::accumulate(x.begin(), x.end(), 0.0) / n;
    double my = std::accumulate(y.begin(), y.end(), 0.0) / n;
    double num = 0.0, den = 0.0;
    for (int i = 0; i < n; ++i) { double dx = x[i]-mx, dy = y[i]-my; num += dx*dy; den += dx*dx; }
    return num/den;
}

// ---------- nice printing helpers ----------
static void print_hr(char c='-', int n=78) { for (int i=0;i<n;++i) std::cout<<c; std::cout<<"\n"; }
static std::string sci(double v, int prec=3) {
    std::ostringstream os; os<<std::scientific<<std::setprecision(prec)<<v; return os.str();
}

int main() {
    const double a = 1.0, T = 1.0, cfl = 0.35;

    // your lists
    std::vector<int> p_list  = {1, 2, 4, 8, 16};
    std::vector<int> ne_list = {16, 32, 64, 128, 256};  // NP (elements)

    // ---------- h-refine for each p: slope log10(L2) vs log10(NP) ----------
    std::cout << "\nDG Advection 1D Convergence \n";
    print_hr();

    for (int p : p_list) {
        std::vector<double> logNP, logL2;
        std::vector<Metrics> metrics(ne_list.size());

        std::cout << "p = " << p << " (np = " << (p+1) << " DOFs/elem)\n";
        std::cout << "  " << std::left
                  << std::setw(8)  << "ne"
                  << std::setw(14) << "h"
                  << std::setw(18) << "L2"
                  << std::setw(18) << "Linf_nodal" << "\n";
        std::cout << "  " << std::setw(8)  << "--"
                  << std::setw(14) << "--"
                  << std::setw(18) << "--"
                  << std::setw(18) << "--" << "\n";

        for (size_t k = 0; k < ne_list.size(); ++k) {
            int ne = ne_list[k];
            Metrics m = run_once(p, ne, a, T, cfl);
            metrics[k] = m;

            double h = 1.0 / ne;
            std::cout << "  " << std::right
                      << std::setw(8)  << ne
                      << std::setw(14) << sci(h, 3)
                      << std::setw(18) << sci(m.L2, 6)
                      << std::setw(18) << sci(m.Linf, 6) << "\n";

            logNP.push_back(std::log10((double)ne));
            logL2.push_back(std::log10(m.L2));
        }

        double slope_np = slope_ls(logNP, logL2);              // expected ≈ -(p+1)
        std::cout << "  Slope log10(L2) vs log10(NP) = "
                  << std::fixed << std::setprecision(3) << slope_np
                  << "   (expected ~ " << -(p+1) << ")\n";
        print_hr();
    }

    // ---------- optional: p-refine summary at fixed ne ----------
    int ne_fixed = 128;
    std::vector<double> P, logL2p;
    std::cout << "p-refine at ne = " << ne_fixed << "\n";
    std::cout << "  " << std::left
              << std::setw(8)  << "p"
              << std::setw(18) << "L2" << "\n";
    std::cout << "  " << std::setw(8)  << "--"
              << std::setw(18) << "--" << "\n";
    for (int p : p_list) {
        Metrics m = run_once(p, ne_fixed, a, T, cfl);
        std::cout << "  " << std::right
                  << std::setw(8)  << p
                  << std::setw(18) << sci(m.L2, 6) << "\n";
        P.push_back((double)p);
        logL2p.push_back(std::log10(m.L2));
    }
    double slope_p = slope_ls(P, logL2p);   // L2 ≈ C * 10^{slope_p * p}
    double factor_per_degree = std::pow(10.0, slope_p);
    std::cout << "  Fit log10(L2) vs p: slope = " << std::fixed << std::setprecision(3)
              << slope_p << "  -> per-degree factor ~ " << factor_per_degree << "\n";

    return 0;
}

// g++ -std=c++20 \
//   basis.cpp quadrature.cpp element.cpp linalg.cpp solver.cpp \
//   examples/main_dg_convergence.cpp \
//   -Iinclude -lopenblas -o main_dg_convergence.exe

// ./ main_dg_convergence.exe
/**
------------------------------------------------------------------------------
p = 1 (np = 2 DOFs/elem)
  ne      h             L2                Linf_nodal
  --      --            --                --
        16     6.250e-02      6.307024e-03      1.765739e-02
        32     3.125e-02      1.387079e-03      3.847870e-03
        64     1.562e-02      3.328178e-04      8.843214e-04
       128     7.812e-03      8.228670e-05      2.110340e-04
       256     3.906e-03      2.051338e-05      5.148444e-05
  Slope log10(L2) vs log10(NP) = -2.060   (expected ~ -2)
------------------------------------------------------------------------------
p = 2 (np = 3 DOFs/elem)
  ne      h             L2                Linf_nodal
  --      --            --                --
        16     6.250e-02      2.592648e-04      9.838899e-04
        32     3.125e-02      3.253231e-05      1.254056e-04
        64     1.562e-02      4.070588e-06      1.574743e-05
       128     7.812e-03      5.089514e-07      1.970608e-06
       256     3.906e-03      6.362294e-08      2.463932e-07
  Slope log10(L2) vs log10(NP) = -2.998   (expected ~ -3)
------------------------------------------------------------------------------
p = 4 (np = 5 DOFs/elem)
  ne      h             L2                Linf_nodal
  --      --            --                --
        16     6.250e-02      1.142009e-07      5.810951e-07
        32     3.125e-02      3.927792e-09      1.860257e-08
        64     1.562e-02      1.261899e-10      6.118629e-10
       128     7.812e-03      3.929631e-12      1.953873e-11
       256     3.906e-03      1.575122e-13      6.458722e-13
  Slope log10(L2) vs log10(NP) = -4.890   (expected ~ -5)
------------------------------------------------------------------------------
p = 8 (np = 9 DOFs/elem)
  ne      h             L2                Linf_nodal
  --      --            --                --
        16     6.250e-02      1.575017e-10      2.227609e-10
        32     3.125e-02      9.868881e-12      1.395668e-11
        64     1.562e-02      6.183329e-13      8.745504e-13
       128     7.812e-03      7.153844e-14      1.103562e-13
       256     3.906e-03      1.206757e-13      1.855183e-13
  Slope log10(L2) vs log10(NP) = -2.781   (expected ~ -9)
------------------------------------------------------------------------------
p = 16 (np = 17 DOFs/elem)
  ne      h             L2                Linf_nodal
  --      --            --                --
        16     6.250e-02      1.112909e-11      1.574305e-11
        32     3.125e-02      6.961005e-13      9.866565e-13
        64     1.562e-02      5.267706e-14      9.203749e-14
       128     7.812e-03      5.804326e-14      1.167955e-13
       256     3.906e-03      1.156871e-13      2.038369e-13
  Slope log10(L2) vs log10(NP) = -1.676   (expected ~ -17)
/