#pragma once
#include <vector>

// Gauss–Legendre on [-1,1] with precision >= p.
// Returns nodes x and weights w using Golub–Welsch via LAPACKE_dstev.
void gauss_legendre(int p, std::vector<double>& x, std::vector<double>& w);
