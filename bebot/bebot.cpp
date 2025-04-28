#include "../include/bebot.h"
#include "../include/bernsteindifferentialmatrix.h"
#include "../include/degelevmatrix.h"
#include "mkl.h"
#include <iostream>
#include <iomanip>

Bebot::Bebot(int N, double T) : N(N), T(T) {}

void Bebot::calculate() {
    // tnodes calculation
    tnodes.resize(N + 1);
    double step = T / N;
    for (int i = 0; i <= N; ++i) {
        tnodes[i] = i * step;
    }



    // w calculation
    w.resize(N + 1, T / (N + 1));



    // dm calculation
    std::vector<double> Dm_temp_flat = BernsteinDifferentiationMatrix(N, T);
    std::vector<double> ElevMatrix_flat = DegElevMatrix(N - 1, N);

    // Initializing size of Dm
    Dm.resize((N + 1) * (N + 1), 0.0);

    // Perform matrix multiplication using MKL
    cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, 
                N + 1, N + 1, N + 1, 
                1.0, 
                &Dm_temp_flat[0], N + 1, 
                &ElevMatrix_flat[0], N + 1, 
                0.0, 
                &Dm[0], N + 1);

}

std::vector<double> Bebot::getNodes() {
    return tnodes;
}

std::vector<double> Bebot::getWeights() {
    return w;
}

const std::vector<double>& Bebot::getDifferentiationMatrix() const {
    return Dm;
}
