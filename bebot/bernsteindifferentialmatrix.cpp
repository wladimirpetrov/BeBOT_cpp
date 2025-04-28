#include "../include/bernsteindifferentialmatrix.h"
#include <vector>

std::vector<double> BernsteinDifferentiationMatrix(int N, double T) {
    int size = (N + 1) * (N + 1);
    std::vector<double> Dm_flat(size, 0.0);

    double factor = static_cast<double>(N) / T; // Explicit cast for precision
    for (int i = 0; i < N; ++i) {
        Dm_flat[i * (N + 1) + i] = -factor;
        if (i + 1 < N + 1) {
            Dm_flat[(i + 1) * (N + 1) + i] = factor;
        }
    }

    return Dm_flat;
}
