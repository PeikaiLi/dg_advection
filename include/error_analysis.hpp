#pragma once

#include <vector>
/*
  L2 error via the mass matrix:
  In each element K we represent the DG function in a nodal Lagrange basis
      u_h(r) = sum_{j=0}^p u_j φ_j(r),  r ∈ [-1,1].
  The continuous L2 inner-product in this finite-dimensional space is represented
  EXACTLY by the (element) mass matrix:
      (M_ref)_{ij} = ∫_{-1}^1 φ_i φ_j dr,   M_el = (h/2) * M_ref.
  Therefore, for two coefficient vectors v,w (i.e., nodal values),
      ∫_K v_h w_h dx  =  v^T M_el w.
  Taking v = w = e (the nodal error/interpolant), we get
      ||e_h||_{L2(K)}^2 = e^T M_el e.
  This matches the DG inner product by construction and is mathematically exact
  for polynomials up to degree 2p (since M_ref was built with Gauss–Legendre).
*/
// L2(error)^2 = sum_e e_e^T * Mel * e_e
double l2_err2_mass(const std::vector<double>& err, int ne, int np, const std::vector<double>& Mel);


// Nodal max-norm: max over stored DOF locations (Chebyshev–GLL points).
// This is not the continuous sup-norm on [0,1], but it’s the standard DG
// "nodal L∞" and is what you can compute directly from the DOFs.
double linf_nodal(const std::vector<double>& err);


// ---------- one run ----------
struct Metrics { double L2; double Linf; };
Metrics run_once(int p, int ne, double a, double T, double cfl, bool fix_time_step = false, double dt_fix = 2*1e-4);