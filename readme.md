# DG Advection (1D Discontinuous Galerkin, Chebyshev nodes, LAPACK/BLAS backend)

This project implements a **Discontinuous Galerkin (DG)** solver for the 1D linear advection equation:

$$
u_t + u_x = 0, \quad x \in [0,1], \quad u(x,0) = u_0(x),
$$

using **Chebyshev nodes** for element discretization and **Gauss–Legendre quadrature** (Golub–Welsch algorithm) for integration.  
Linear algebra operations are performed through **LAPACK/BLAS** (e.g. `dgemm`, `dgetrf`, `dgetrs`, `dstev`).

---

## Features
- Discontinuous Galerkin method with polynomial basis of order `p`
- Chebyshev nodes mapping for local interpolation
- Gauss–Legendre quadrature (computed via Golub–Welsch)
- Mass and stiffness matrix assembly per element
- Explicit Runge–Kutta 4 (RK4) time stepping
- Periodic boundary conditions, upwind flux
- Error computation (∞-norm and L²-norm)
- Convergence rate study

---

## Project Structure
All sources are kept in a flat layout:

```text
dg_advection/
├── Makefile             # Build configuration (GNU Make)
├── main.cpp             # Entry point, runs demo + convergence study
├── basis.hpp / .cpp     # Legendre basis evaluation and derivatives
├── quadrature.hpp / .cpp# Gauss–Legendre quadrature (Golub–Welsch)
├── element.hpp / .cpp   # Element-level matrices (M, C) and interpolation
├── solver.hpp / .cpp    # DG RHS assembly, Runge–Kutta time stepping
├── analysis.hpp / .cpp  # Error norms, convergence slope computation
├── linalg.hpp / .cpp    # Lightweight matrix wrapper, BLAS/LAPACK interface
```

This keeps everything simple in one directory.  
If the project grows, it can later be reorganized into `include/` + `src/`.

---

## Dependencies
- **C++17** compiler (I am using MSVC with a GNU Make environment)
- **GNU Make**
- **BLAS/LAPACK** implementation (choose one):
  - [OpenBLAS](https://www.openblas.net/) (recommended, cross-platform)
  - Intel oneMKL (alternative on Windows/Linux)

On Windows, the easiest setup is [MSYS2 with MinGW-w64](https://www.msys2.org/).  

---

## Build Instructions (with Make)

1. Clone or download the project, and enter the folder:

   ```bash
   cd dg_advection
   ```

2. Build the executable:

   ```bash
   make
   ```

   This compiles all `.cpp` files and links against BLAS/LAPACK.

3. Clean build artifacts:

   ```bash
   make clean
   ```

---

## Run
Execute the solver:

```bash
./dg_advection.exe
```

Expected output:
- Error norms (inf-norm, two-norm) for a test case
- Convergence rates for different polynomial orders

---

## Next Steps
- Add Python scripts (matplotlib) to plot/export results

---

## License
MIT License — free to use for research and education.
