#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <numbers>
#include <numeric>
#include "solver.hpp"
#include "error_analysis.hpp"


// ---------- nice printing helpers ----------
static void print_hr(char c='-', int n=78) { for (int i=0;i<n;++i) std::cout<<c; std::cout<<"\n"; }
static std::string sci(double v, int prec=3) {
    std::ostringstream os; os<<std::scientific<<std::setprecision(prec)<<v; return os.str();
}

int main() {
    const double a = 1.0, T = 1.0, cfl = 0.35;

    std::vector<int> p_list  = {1, 2, 4, 8, 16};
    std::vector<int> ne_list = {16, 32, 64, 128, 256};  // target total DOFs

    // errors[p_index][ne_index]
    std::vector<std::vector<double>> errors(p_list.size(),
                                            std::vector<double>(ne_list.size(), 0.0));

    std::cout << "\nDG Advection 1D Convergence \n";
    print_hr();

    for (size_t ip = 0; ip < p_list.size(); ++ip) {
        int p = p_list[ip];
        std::vector<double> logNP, logL2;

        std::cout << "p = " << p << " (np = " << (p+1) << " DOFs/elem)\n";
        std::cout << "        " << std::left
                  << std::setw(8)  << "ne"
                  << std::setw(14) << "h"
                  << std::setw(18) << "L2"
                  << std::setw(18) << "Linf_nodal" << "\n";
        std::cout << "        " << std::setw(8)  << "--"
                  << std::setw(14) << "--"
                  << std::setw(18) << "--"
                  << std::setw(18) << "--" << "\n";

        for (size_t k = 0; k < ne_list.size(); ++k) {
            int ne = ne_list[k] / (p);  // number of elements
            // Metrics m = run_once_DGAdvection1D(p, ne, a, T, cfl, true, 2*1e-4);
            // true for fixed time step

            // automatic time-stepping
            Metrics m = run_once_DGAdvection1D(p, ne, a, T, cfl, false, 2*1e-4);

            double h = 1.0 / ne;
            std::cout << "  " << std::right
                      << std::setw(8)  << ne
                      << std::setw(14) << sci(h, 3)
                      << std::setw(18) << sci(m.L2, 6)
                      << std::setw(18) << sci(m.Linf, 6) << "\n";

            logNP.push_back(std::log10((double)ne));
            logL2.push_back(std::log10(m.L2));

            // out into errors matrix
            errors[ip][k] = m.L2;
        }

        double slope_np = slope_max_adjacent(logNP, logL2);
        std::cout << "  Slope log10(L2) vs log10(NP) = "
                  << std::fixed << std::setprecision(3) << slope_np
                  << "   (expected ~ " << -(p+1) << ")\n";
        print_hr();
    }

    // ---------- print 5×5 L2 error table ----------
    std::cout << "\n5x5 L2 error table (rows=p, cols=ne):\n";
    for (size_t ip = 0; ip < p_list.size(); ++ip) {
        for (size_t k = 0; k < ne_list.size(); ++k) {
            std::cout << std::scientific << std::setprecision(6)
                      << errors[ip][k] << "  ";
        }
        std::cout << "\n";
    }
}

// g++ -std=c++20   basis.cpp quadrature.cpp element.cpp linalg.cpp solver.cpp solver_advdiff.cpp error_analysis.cpp  examples/main_dg_convergence.cpp   -Iinclude -lopenblas -o main_dg_convergence.exe
// ./main_dg_convergence.exe

// fixed time step
// DG Advection 1D Convergence 
// ------------------------------------------------------------------------------
// p = 1 (np = 2 DOFs/elem)
//         ne      h             L2                Linf_nodal        
//         --      --            --                --
//         16     6.250e-02      4.723952e-02      1.889903e-01
//         32     3.125e-02      1.196027e-02      5.156706e-02
//         64     1.562e-02      1.941683e-03      9.868625e-03
//        128     7.812e-03      2.973049e-04      1.793201e-03
//        256     3.906e-03      5.403365e-05      3.529185e-04
//   Slope log10(L2) vs log10(NP) = -2.707   (expected ~ -2)
// ------------------------------------------------------------------------------
// p = 2 (np = 3 DOFs/elem)
//         ne      h             L2                Linf_nodal
//         --      --            --                --
//          8     1.250e-01      2.911725e-02      1.199747e-01
//         16     6.250e-02      3.825109e-03      1.476169e-02
//         32     3.125e-02      2.849124e-04      1.779697e-03
//         64     1.562e-02      3.197235e-05      2.482672e-04
//        128     7.812e-03      3.975629e-06      3.100645e-05
//   Slope log10(L2) vs log10(NP) = -3.747   (expected ~ -3)
// ------------------------------------------------------------------------------
// p = 4 (np = 5 DOFs/elem)
//         ne      h             L2                Linf_nodal
//         --      --            --                --
//          4     2.500e-01      2.255224e-02      5.897820e-02
//          8     1.250e-01      8.628952e-04      6.054906e-03
//         16     6.250e-02      1.823028e-05      2.160166e-04
//         32     3.125e-02      5.725496e-07      6.488529e-06
//         64     1.562e-02      1.940513e-08      1.962599e-07
//   Slope log10(L2) vs log10(NP) = -5.565   (expected ~ -5)
// ------------------------------------------------------------------------------
// p = 8 (np = 9 DOFs/elem)
//         ne      h             L2                Linf_nodal
//         --      --            --                --
//          2     5.000e-01      1.467562e-02      4.667126e-02
//          4     2.500e-01      1.406573e-04      6.507371e-04
//          8     1.250e-01      2.221917e-07      2.302978e-06
//         16     6.250e-02      8.845158e-10      9.746696e-09
//         32     3.125e-02      1.458055e-11      5.671696e-11
//   Slope log10(L2) vs log10(NP) = -9.306   (expected ~ -9)
// ------------------------------------------------------------------------------
// p = 16 (np = 17 DOFs/elem)
//         ne      h             L2                Linf_nodal
//         --      --            --                --
//          1     1.000e+00      1.331652e-02      2.580537e-02
//          2     5.000e-01      1.747688e-05      4.473550e-05
//          4     2.500e-01      6.566400e-10      6.848505e-09
//          8     1.250e-01      1.451289e-11      4.297984e-11
//         16     6.250e-02      1.451113e-11      4.361489e-11
//   Slope log10(L2) vs log10(NP) = -14.700   (expected ~ -17)
// ------------------------------------------------------------------------------

// 5x5 L2 error table (rows=p, cols=ne):
// 4.723952e-02  1.196027e-02  1.941683e-03  2.973049e-04  5.403365e-05
// 2.911725e-02  3.825109e-03  2.849124e-04  3.197235e-05  3.975629e-06
// 2.255224e-02  8.628952e-04  1.823028e-05  5.725496e-07  1.940513e-08
// 1.467562e-02  1.406573e-04  2.221917e-07  8.845158e-10  1.458055e-11
// 1.331652e-02  1.747688e-05  6.566400e-10  1.451289e-11  1.451113e-11