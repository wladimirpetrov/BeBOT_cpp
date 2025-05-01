#include "../include/piecewisebernsteinproduct.h"
#include "../include/bernsteinproduct.h"
#include <vector>
#include <cmath>
#include <algorithm>

std::vector<double> PiecewiseBernsteinProduct(
    const std::vector<double>& A,
    const std::vector<double>& B,
    int K,
    int N)
{
    // each segment is length (N+1), and BernsteinProduct produces length (2N+1)
    const int rows    = N + 1;
    const int outSize = 2 * N + 1;

    std::vector<double> Cp; 
    Cp.resize(static_cast<size_t>(K) * outSize);

    // temporary buffers to hold one segment at a time
    std::vector<double> Ai(rows), Bi(rows);
    std::vector<double> Ci;

    for (int seg = 0; seg < K; ++seg) {
        // pointer to where this segment's output goes
        double* Cp_ptr = Cp.data() + seg * outSize;

        // copy & sanitize this segment
        for (int j = 0; j < rows; ++j) {
            double a = A[ seg * rows + j ];
            double b = B[ seg * rows + j ];
            Ai[j] = std::isnan(a) ? 0.0 : a;
            Bi[j] = std::isnan(b) ? 0.0 : b;
        }

        // compute product
        // your existing interface returns a new vector:
        Ci = BernsteinProduct(Ai, Bi);

        // just in case, guard against a malformed implementation
        if ((int)Ci.size() != outSize) {
            // fill with zeros instead of returning NaNs
            std::fill_n(Cp_ptr, outSize, 0.0);
        } else {
            // copy the result into the flat output buffer
            std::copy(Ci.begin(), Ci.end(), Cp_ptr);
        }
    }

    return Cp;
}
