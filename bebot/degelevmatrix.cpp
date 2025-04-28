#include "../include/degelevmatrix.h"
#include "../include/nchoosek_mod.h"
#include <vector>

std::vector<double> DegElevMatrix(int N, int M) {
    int r = M - N;
    std::vector<std::vector<double>> E(M + 1, std::vector<double>(N + 1, 0.0));

    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= r; ++j) {
            double numerator = nchoosek_mod(N, i) * nchoosek_mod(r, j);
            double denominator = nchoosek_mod(M, i + j);
            E[i + j][i] = numerator / denominator;
        }
    }

    std::vector<double> E_flattened((M + 1) * (N + 1), 0.0);
    for (size_t i = 0; i < E.size(); ++i) {
        for (size_t j = 0; j < E[i].size(); ++j) {
            E_flattened[j * (M + 1) + i] = E[i][j];
        }
    }

    return E_flattened;
}
