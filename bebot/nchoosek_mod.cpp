#include "../include/nchoosek_mod.h"

// Stable computation of binomial coefficient
double nchoosek_mod(int N, int k) {
    if (k < 0 || k > N) return 0.0;
    if (k == 0 || k == N) return 1.0;

    double log_binom = 0.0;
    for (int i = 1; i <= k; ++i) {
        log_binom += std::log(N - (k - i)) - std::log(i);
    }
    return std::exp(log_binom);
}
