#include "../include/piecewisedegelevmatrix.h"
#include "../include/degelevmatrix.h"
#include <vector>
#include <mkl.h>   // MKL CBLAS header

std::vector<double> PiecewiseDegElevMatrix(
    const std::vector<double>& Cp, // size (N+1)*K, column-major
    int K,
    int N,
    int M
) {
    // 1) build the (M+1)x(N+1) elevation matrix once
    std::vector<double> E_flat = DegElevMatrix(N, M);
    
    // 2) prepare output (M+1)xK
    std::vector<double> Cp_elev((M + 1) * K, 0.0);

    // 3) C ← 1.0·E_flat·Cp + 0.0·C
    cblas_dgemm(
        CblasColMajor,      // our data is column-major
        CblasNoTrans,       // don’t transpose E_flat
        CblasNoTrans,       // don’t transpose Cp
        M+1,                // number of rows of E_flat and Cp_elev
        K,                  // number of cols of Cp and Cp_elev
        N+1,                // shared dimension
        1.0,                // α
        E_flat.data(),      // A
        M+1,                // leading dim of A (rows of E_flat)
        Cp.data(),          // B
        N+1,                // leading dim of B (rows of Cp)
        0.0,                // β
        Cp_elev.data(),     // C
        M+1                 // leading dim of C
    );

    return Cp_elev;
}
