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

            // true for fixed time step
            Metrics m = run_once_DGAdvection1D_Diff(p, ne, 1e-3, a, T, cfl, true, 6e-5);


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

// g++ -std=c++20   basis.cpp quadrature.cpp element.cpp linalg.cpp solver.cpp solver_advdiff.cpp error_analysis.cpp  examples/main_dg_conv_diff_best.cpp   -Iinclude -lopenblas -o main_dg_conv_diff_best.exe
// ./main_dg_conv_diff_best.exe
/*DG Advection 1D Convergence 
------------------------------------------------------------------------------
p = 1 (np = 2 DOFs/elem)
        ne      h             L2                Linf_nodal        
        --      --            --                --
        16     6.250e-02      2.986319e-02      1.113763e-01
        32     3.125e-02      6.286708e-03      2.488510e-02
        64     1.562e-02      9.029511e-04      3.772525e-03
       128     7.812e-03      1.133650e-04      4.238191e-04
       256     3.906e-03      1.540331e-05      7.527928e-05
  Slope log10(L2) vs log10(NP) = -2.994   (expected ~ -2)
------------------------------------------------------------------------------
p = 2 (np = 3 DOFs/elem)
        ne      h             L2                Linf_nodal        
        --      --            --                --
         8     1.250e-01      1.644297e-02      5.915699e-02
        16     6.250e-02      1.546910e-03      6.350989e-03
        32     3.125e-02      1.182970e-04      7.572652e-04
        64     1.562e-02      1.075520e-05      7.631784e-05
       128     7.812e-03      8.702816e-07      5.759978e-06
  Slope log10(L2) vs log10(NP) = -3.709   (expected ~ -3)
------------------------------------------------------------------------------
p = 4 (np = 5 DOFs/elem)
        ne      h             L2                Linf_nodal
        --      --            --                --
         4     2.500e-01      1.059208e-02      2.379032e-02
         8     1.250e-01      2.946851e-04      2.165101e-03
        16     6.250e-02      6.579047e-06      4.566865e-05
        32     3.125e-02      1.845055e-07      1.096440e-06
        64     1.562e-02      4.901948e-09      1.917845e-08
  Slope log10(L2) vs log10(NP) = -5.485   (expected ~ -5)
------------------------------------------------------------------------------
p = 8 (np = 9 DOFs/elem)
        ne      h             L2                Linf_nodal
        --      --            --                --
         2     5.000e-01      5.183391e-03      1.433837e-02
         4     2.500e-01      2.519780e-05      9.282344e-05
         8     1.250e-01      6.223130e-08      3.561629e-07
        16     6.250e-02      1.197279e-10      6.128786e-10
        32     3.125e-02      6.993029e-13      2.331892e-12
  Slope log10(L2) vs log10(NP) = -9.022   (expected ~ -9)
------------------------------------------------------------------------------
p = 16 (np = 17 DOFs/elem)
        ne      h             L2                Linf_nodal
        --      --            --                --
         1     1.000e+00      4.962711e-03      1.072450e-02
         2     5.000e-01      1.213491e-06      2.648716e-06
         4     2.500e-01      2.396085e-11      9.818995e-11
         8     1.250e-01      6.716694e-13      2.331963e-12
        16     6.250e-02      6.716568e-13      2.331964e-12
  Slope log10(L2) vs log10(NP) = -15.628   (expected ~ -17)
------------------------------------------------------------------------------

5x5 L2 error table (rows=p, cols=ne):
2.986319e-02  6.286708e-03  9.029511e-04  1.133650e-04  1.540331e-05
1.644297e-02  1.546910e-03  1.182970e-04  1.075520e-05  8.702816e-07
1.059208e-02  2.946851e-04  6.579047e-06  1.845055e-07  4.901948e-09
5.183391e-03  2.519780e-05  6.223130e-08  1.197279e-10  6.993029e-13
4.962711e-03  1.213491e-06  2.396085e-11  6.716694e-13  6.716568e-13
 */
