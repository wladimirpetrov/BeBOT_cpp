#ifndef PIEWISEDEGELEVMATRIX_H
#define PIEWISEDEGELEVMATRIX_H

#include <vector>

// Piecewise degree elevation for a K-segment Bernstein curve.
// Cp: flattened control points of length (N+1)*K, laid out segment-by-segment.
// K : number of segments
// N : original polynomial degree
// M : desired elevated degree
// Returns a flattened vector of length (M+1)*K.
std::vector<double> PiecewiseDegElevMatrix(
    const std::vector<double>& Cp,
    int K,
    int N,
    int M
);

#endif

