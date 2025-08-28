#pragma once
#include <vector>

// Build reference operators on [-1,1] for given polynomial degree p.
// Outputs:
//   s      : Chebyshev–G-L nodes (size p+1, descending)
//   gx, gw : Gauss–Legendre nodes/weights (size p+1, exact for degree 2p)
//   phi_Q  : values of nodal basis at gx,   size (ng x (p+1)) row-major
//   dphi_Q : derivatives of nodal basis at gx,   size (ng x (p+1)) row-major
//   Mref   : reference mass matrix  (p+1 x p+1), row-major
//   Cref   : reference convection   (p+1 x p+1), row-major
void build_reference_ops(int p,
                         std::vector<double>& s,
                         std::vector<double>& gx, std::vector<double>& gw,
                         std::vector<double>& phi_Q, std::vector<double>& dphi_Q,
                         std::vector<double>& Mref,  std::vector<double>& Cref);
