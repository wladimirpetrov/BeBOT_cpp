#include "state_space.h"
#include <vector>

std::vector<std::vector<double>> get_A_matrix(double speed) {
    // Define A matrices for different speeds. For simplicity, only two are shown.
    std::vector<std::vector<double>> A1 = {
        {-9.61335633071093e-13, 4.99999999991687, 1.00000000005217, 1.12249207912694e-10},
        {1.45212601715332e-12, -5.13981176158834e-14, -1.17800813245335e-10, 1.00000000101190},
        {0.000840272043916469, -3.04085570568224e-07, -0.101668288331192, -2.70007795164826},
        {4.37929332288747e-06, -0.0149473707639360, -0.00572337288021148, -0.244530705224609}
    };
    std::vector<std::vector<double>> A2 = {
        // Define another A matrix for a different speed.
    };

    if (speed < 10.0) {
        return A1;
    } else {
        return A2;
    }
}

std::vector<std::vector<double>> get_B_matrix(double speed) {
    // Define B matrices for different speeds. For simplicity, only two are shown.
    std::vector<std::vector<double>> B1 = {
        {-1.69674752574204e-12, -1.10484134201392e-17, 2.75706811981441e-12},
        {-7.22760359972355e-13, -6.12418309548146e-18, -4.20492416927654e-12},
        {-0.00261112713870843, -2.21991902099372e-06, 0.00164151357821695},
        {0.000330132364452263, 2.40305021375118e-16, 7.27066721427474e-05}
    };
    std::vector<std::vector<double>> B2 = {
        // Define another B matrix for a different speed.
    };

    if (speed < 10.0) {
        return B1;
    } else {
        return B2;
    }
}
