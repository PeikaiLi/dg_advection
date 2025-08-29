#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <numbers>
#include "solver.hpp"
#include "error_analysis.hpp"

int main() {
    const int    p   = 16;
    const int    Ne  = 16;     // number of elements
    const double a   = 1.0;   // advection speed
    const double cfl = 0.35;
    const double T   = 1.0;   // one period for sin(2πx) with a=1

    std::cout << "start...\n";
    DGAdvection1D dg(p, Ne, a);
    std::cout << "finished dg constructor\n";

    // solve the PDE
    dg.solve_pde(cfl, T);
    return 0;
}

// g++ -std=c++20 \
//   basis.cpp quadrature.cpp element.cpp linalg.cpp solver.cpp error_analysis.cpp \
//   examples/check_advection.cpp \
//   -Iinclude -lopenblas -o check_advection.exe
// ./check_advection.exe

// start...
// finished dg constructor
// p=16, ne=16, dt=0.000662691, steps=1509
//   L2 error squared (mass-matrix): 1.238566e-22
//   L2 error:                       1.112909e-11
//   Linf error (nodal):             1.574305e-11