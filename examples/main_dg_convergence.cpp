#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <numbers>
#include <numeric>
#include "solver.hpp"
#include "error_analysis.hpp"



// // ---------- least-squares slope y = m x + b ----------
// static double slope_ls(const std::vector<double>& x, const std::vector<double>& y) {
//     const int n = (int)x.size();
//     double mx = std::accumulate(x.begin(), x.end(), 0.0) / n;
//     double my = std::accumulate(y.begin(), y.end(), 0.0) / n;
//     double num = 0.0, den = 0.0;
//     for (int i = 0; i < n; ++i) { double dx = x[i]-mx, dy = y[i]-my; num += dx*dy; den += dx*dx; }
//     return num/den;
// }

// slope between the two neighboring points with maximum slope (by absolute value)
static double slope_max_adjacent(const std::vector<double>& x,
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
        std::cout << "         " << std::left
                  << std::setw(8)  << "ne"
                  << std::setw(14) << "h"
                  << std::setw(18) << "L2"
                  << std::setw(18) << "Linf_nodal" << "\n";
        std::cout << "         " << std::setw(8)  << "--"
                  << std::setw(14) << "--"
                  << std::setw(18) << "--"
                  << std::setw(18) << "--" << "\n";

        for (size_t k = 0; k < ne_list.size(); ++k) {
            int ne = ne_list[k] / (p);  // number of elements
            Metrics m = run_once(p, ne, a, T, cfl, false, 2*1e-4);

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

// g++ -std=c++20   basis.cpp quadrature.cpp element.cpp linalg.cpp solver.cpp error_analysis.cpp  examples/main_dg_convergence.cpp   -Iinclude -lopenblas -o main_dg_convergence.exe
// ./main_dg_convergence.exe

// fixed time step
// Peikai_Li@LPK UCRT64 /d/UCB/course/Math228B/pj/PJ/PS7DG_cpp/dg_cpp_demo/dg_advection
// $ ./main_dg_convergence.exe

// DG Advection 1D Convergence
// ------------------------------------------------------------------------------
// p = 1 (np = 2 DOFs/elem)
//            ne      h             L2                Linf_nodal
//            --      --            --                --
//         16     6.250e-02      6.307030e-03      1.765739e-02
//         32     3.125e-02      1.387079e-03      3.847870e-03
//         64     1.562e-02      3.328178e-04      8.843214e-04
//        128     7.812e-03      8.228670e-05      2.110340e-04
//        256     3.906e-03      2.051338e-05      5.148444e-05
//   Slope log10(L2) vs log10(NP) = -2.185   (expected ~ -2)
// ------------------------------------------------------------------------------
// p = 2 (np = 3 DOFs/elem)
//            ne      h             L2                Linf_nodal
//            --      --            --                --
//          8     1.250e-01      2.047506e-03      7.208359e-03
//         16     6.250e-02      2.592648e-04      9.839176e-04
//         32     3.125e-02      3.253231e-05      1.254074e-04
//         64     1.562e-02      4.070588e-06      1.574755e-05
//        128     7.812e-03      5.089514e-07      1.970615e-06
//   Slope log10(L2) vs log10(NP) = -3.000   (expected ~ -3)
// ------------------------------------------------------------------------------
// p = 4 (np = 5 DOFs/elem)
//            ne      h             L2                Linf_nodal
//            --      --            --                --
//          4     2.500e-01      1.113473e-04      5.154091e-04
//          8     1.250e-01      3.448187e-06      1.844945e-05
//         16     6.250e-02      1.142215e-07      5.783492e-07
//         32     3.125e-02      3.927272e-09      1.842645e-08
//         64     1.562e-02      1.259546e-10      6.009691e-10
//   Slope log10(L2) vs log10(NP) = -5.013   (expected ~ -5)
// ------------------------------------------------------------------------------
// p = 8 (np = 9 DOFs/elem)
//            ne      h             L2                Linf_nodal
//            --      --            --                --
//          2     5.000e-01      6.003355e-08      3.106835e-07
//          4     2.500e-01      9.351614e-10      4.740444e-09
//          8     1.250e-01      1.839369e-12      1.188472e-11
//         16     6.250e-02      9.244620e-14      1.548782e-13
//         32     3.125e-02      9.354770e-14      1.326717e-13
//   Slope log10(L2) vs log10(NP) = -8.990   (expected ~ -9)
// ------------------------------------------------------------------------------
// p = 16 (np = 17 DOFs/elem)
//            ne      h             L2                Linf_nodal
//            --      --            --                --
//          1     1.000e+00      1.048974e-11      6.979242e-11
//          2     5.000e-01      9.320860e-14      1.332753e-13
//          4     2.500e-01      9.308788e-14      1.414025e-13
//          8     1.250e-01      9.308373e-14      1.391031e-13
//         16     6.250e-02      9.300598e-14      1.358262e-13
//   Slope log10(L2) vs log10(NP) = -6.814   (expected ~ -17)
// ------------------------------------------------------------------------------

// 5x5 L2 error table (rows=p, cols=ne):
// 6.307030e-03  1.387079e-03  3.328178e-04  8.228670e-05  2.051338e-05
// 2.047506e-03  2.592648e-04  3.253231e-05  4.070588e-06  5.089514e-07
// 1.113473e-04  3.448187e-06  1.142215e-07  3.927272e-09  1.259546e-10
// 6.003355e-08  9.351614e-10  1.839369e-12  9.244620e-14  9.354770e-14
// 1.048974e-11  9.320860e-14  9.308788e-14  9.308373e-14  9.300598e-14



/* automatic time step with CFL
DG Advection 1D Convergence 
------------------------------------------------------------------------------
p = 1 (np = 2 DOFs/elem)
           ne      h             L2                Linf_nodal        
           --      --            --                --
        16     6.250e-02      6.307024e-03      1.765739e-02
        32     3.125e-02      1.387079e-03      3.847870e-03
        64     1.562e-02      3.328178e-04      8.843214e-04
       128     7.812e-03      8.228670e-05      2.110340e-04
       256     3.906e-03      2.051338e-05      5.148444e-05
  Slope log10(L2) vs log10(NP) = -2.185   (expected ~ -2)
------------------------------------------------------------------------------
p = 2 (np = 3 DOFs/elem)
           ne      h             L2                Linf_nodal
           --      --            --                --
         8     1.250e-01      2.047503e-03      7.207891e-03
        16     6.250e-02      2.592648e-04      9.838899e-04
        32     3.125e-02      3.253231e-05      1.254056e-04
        64     1.562e-02      4.070588e-06      1.574743e-05
       128     7.812e-03      5.089514e-07      1.970608e-06
  Slope log10(L2) vs log10(NP) = -3.000   (expected ~ -3)
------------------------------------------------------------------------------
p = 4 (np = 5 DOFs/elem)
           ne      h             L2                Linf_nodal
           --      --            --                --
         4     2.500e-01      1.113211e-04      5.153484e-04
         8     1.250e-01      3.447498e-06      1.848335e-05
        16     6.250e-02      1.142009e-07      5.810951e-07
        32     3.125e-02      3.927792e-09      1.860257e-08
        64     1.562e-02      1.261899e-10      6.118629e-10
  Slope log10(L2) vs log10(NP) = -5.013   (expected ~ -5)
------------------------------------------------------------------------------
p = 8 (np = 9 DOFs/elem)
           ne      h             L2                Linf_nodal
           --      --            --                --
         2     5.000e-01      6.284474e-07      9.178877e-07
         4     2.500e-01      3.991898e-08      6.080590e-08
         8     1.250e-01      2.520019e-09      3.575126e-09
        16     6.250e-02      1.575017e-10      2.227609e-10
        32     3.125e-02      9.868881e-12      1.395668e-11
  Slope log10(L2) vs log10(NP) = -4.000   (expected ~ -9)
------------------------------------------------------------------------------
p = 16 (np = 17 DOFs/elem)
           ne      h             L2                Linf_nodal
           --      --            --                --
         1     1.000e+00      7.084184e-07      1.001844e-06
         2     5.000e-01      4.522214e-08      6.395357e-08
         4     2.500e-01      2.826408e-09      3.997140e-09
         8     1.250e-01      1.775892e-10      2.511490e-10
        16     6.250e-02      1.112909e-11      1.574305e-11
  Slope log10(L2) vs log10(NP) = -4.000   (expected ~ -17)
------------------------------------------------------------------------------

5x5 L2 error table (rows=p, cols=ne):
6.307024e-03  1.387079e-03  3.328178e-04  8.228670e-05  2.051338e-05
2.047503e-03  2.592648e-04  3.253231e-05  4.070588e-06  5.089514e-07
1.113211e-04  3.447498e-06  1.142009e-07  3.927792e-09  1.261899e-10
6.284474e-07  3.991898e-08  2.520019e-09  1.575017e-10  9.868881e-12
7.084184e-07  4.522214e-08  2.826408e-09  1.775892e-10  1.112909e-11
*/