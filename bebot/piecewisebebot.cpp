#include "../include/piecewisebebot.h"
#include "../include/bernsteindifferentialmatrix.h"
#include "../include/degelevmatrix.h"
#include <vector>
#include <cmath>
#include <limits>

PiecewiseBeBOT::PiecewiseBeBOT(int N, const std::vector<double>& tknots)
    : N(N), originalTknots(tknots) {
    calculate();
}

void PiecewiseBeBOT::transformTknots() {
    transformedTknots.clear();
    for (size_t i = 0; i + 1 < originalTknots.size(); ++i) {
        transformedTknots.emplace_back(
            std::vector<double>{ originalTknots[i], originalTknots[i+1] } );
    }
}

void PiecewiseBeBOT::calculate() {
    transformTknots();
    const double T = originalTknots.back() - originalTknots.front();
    const int M = static_cast<int>(transformedTknots.size());
    const int dim = N + 1;

    // allocate outputs
    tnodes.resize(dim * M);
    w.assign(dim * M, T / (dim * M));
    Dm_flat.assign(dim * dim * M, 0.0);

    for (int seg = 0; seg < M; ++seg) {
        double t0 = transformedTknots[seg][0];
        double t1 = transformedTknots[seg][1];
        double segLen = t1 - t0;
        double eps = std::numeric_limits<double>::epsilon() * segLen;

        // compute local nodes with MATLAB-style eps adjustment
        double a = (seg == 0 ? t0 : t0 + eps);
        double b = (seg == M-1 ? t1 : t1 - eps);
        for (int i = 0; i < dim; ++i) {
            tnodes[seg * dim + i] = a + (b - a) * (static_cast<double>(i) / N);
        }

        // get differentiation and degree-elevation matrices
        std::vector<double> Dm_temp = BernsteinDifferentiationMatrix(N, segLen);  // size dim×dim
        std::vector<double> Elev    = DegElevMatrix(N-1, N);                     // size dim×dim

        // compute block = Dm_temp * Elev into blockDm
        std::vector<double> blockDm(dim * dim, 0.0);
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dim; ++k) {
                    sum += Dm_temp[i*dim + k] * Elev[k*dim + j];
                }
                blockDm[i*dim + j] = sum;
            }
        }

        // copy blockDm into block diagonal of Dm_flat
        size_t base = static_cast<size_t>(seg) * dim * dim;
        for (size_t k = 0; k < static_cast<size_t>(dim*dim); ++k) {
            Dm_flat[base + k] = blockDm[k];
        }
    }
}

std::vector<double> PiecewiseBeBOT::getNodes() const {
    return tnodes;
}

std::vector<double> PiecewiseBeBOT::getWeights() const {
    return w;
}

std::vector<std::vector<double>> PiecewiseBeBOT::getDifferentiationMatrix() const {
    const int M = static_cast<int>(transformedTknots.size());
    const int dim = N + 1;
    const int total = dim * M;
    std::vector<std::vector<double>> Dm2D(total, std::vector<double>(total, 0.0));
    for (int seg = 0; seg < M; ++seg) {
        size_t base = static_cast<size_t>(seg) * dim * dim;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                Dm2D[seg*dim + i][seg*dim + j] = Dm_flat[base + i*dim + j];
            }
        }
    }
    return Dm2D;
}

const std::vector<double>& PiecewiseBeBOT::getDifferentiationMatrixFlat() const {
    return Dm_flat;
}
