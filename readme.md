# DG Advection – 1D Discontinuous Galerkin Solver

A C++ implementation of high–order **Discontinuous Galerkin (DG)** methods for 1D **advection** and **advection–diffusion** equations with periodic BCs. 

## Features
- Polynomial DG(p) discretization, Chebyshev nodes
- Gauss–Legendre quadrature (Golub–Welsch)
- Upwind flux (advection), LDG scheme (diffusion)
- Explicit RK4 time integration
- BLAS/LAPACK (OpenBLAS) for matrix ops
- Error analysis: L² norm, nodal L∞ norm, convergence rates

## Structure
```
basis.*        # Basis functions
quadrature.*   # Gauss–Legendre quadrature
element.*      # Element operators
linalg.*       # BLAS/LAPACK wrappers
solver.*       # Advection solver
solver_advdiff.* # Advection–diffusion (LDG)
error_analysis.* # Error metrics
examples/      # Driver programs
```
Executables such as `main_dg_convergence`, `main_dgconv_diff`, `main_dg_conv_diff_best` run tests.

## Build
### Windows (MSYS2 UCRT64)
```bash
pacman -Syu
pacman -S mingw-w64-ucrt-x86_64-toolchain mingw-w64-ucrt-x86_64-openblas mingw-w64-ucrt-x86_64-lapack 
```
### Linux/macOS
```bash
g++ -std=c++20 *.cpp -Iinclude -lopenblas -o dg_advection
```

## Usage
```bash
./main_dg_convergence
./main_dgconv_diff
./main_dg_conv_diff_best
```

## Notes
- CFL condition: Δt ≲ h / ((2p+1)|a|)  is not good for diffusion equations
- Errors decrease at ~p+1 order for smooth solutions

## License
MIT
