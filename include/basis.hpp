#pragma once
#include <vector>

// Return Chebyshev–Gauss–Lobatto nodes 
// s_j = cos(j*pi/p), 
// j=p..0 (descending so that s_j goes from -1 to 1)
std::vector<double> chebyshev_nodes(int p);

// Evaluate Legendre polynomials P_0..P_p 
// and their derivatives at points x.
// Y has size (nx*(p+1)), row-major: Y[i*(p+1)+k] = P_k(x_i)
// dY same layout for P'_k(x_i)
void legendre_poly(const std::vector<double>& x, 
                   const int p,
                   std::vector<double>& Y,
                   std::vector<double>& dY);
