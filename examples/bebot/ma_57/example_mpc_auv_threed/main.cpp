// main.cpp
#include "state_space_matrices.h"
#include <iostream>
#include <vector>

void extractAndPrintMatrixElements(const std::vector<double>& A, const std::vector<double>& B) {
    // Assuming A is 4x4
    double a11 = A[0], a12 = A[1], a13 = A[2], a14 = A[3];
    double a21 = A[4], a22 = A[5], a23 = A[6], a24 = A[7];
    double a31 = A[8], a32 = A[9], a33 = A[10], a34 = A[11];
    double a41 = A[12], a42 = A[13], a43 = A[14], a44 = A[15];

    std::cout << "Matrix A elements:" << std::endl;
    std::cout << "a11 = " << a11 << ", a12 = " << a12 << ", a13 = " << a13 << ", a14 = " << a14 << std::endl;
    std::cout << "a21 = " << a21 << ", a22 = " << a22 << ", a23 = " << a23 << ", a24 = " << a24 << std::endl;
    std::cout << "a31 = " << a31 << ", a32 = " << a32 << ", a33 = " << a33 << ", a34 = " << a34 << std::endl;
    std::cout << "a41 = " << a41 << ", a42 = " << a42 << ", a43 = " << a43 << ", a44 = " << a44 << std::endl;

    // Assuming B is 4x3
    double b11 = B[0], b12 = B[1], b13 = B[2];
    double b21 = B[3], b22 = B[4], b23 = B[5];
    double b31 = B[6], b32 = B[7], b33 = B[8];
    double b41 = B[9], b42 = B[10], b43 = B[11];

    std::cout << "Matrix B elements:" << std::endl;
    std::cout << "b11 = " << b11 << ", b12 = " << b12 << ", b13 = " << b13 << std::endl;
    std::cout << "b21 = " << b21 << ", b22 = " << b22 << ", b23 = " << b23 << std::endl;
    std::cout << "b31 = " << b31 << ", b32 = " << b32 << ", b33 = " << b33 << std::endl;
    std::cout << "b41 = " << b41 << ", b42 = " << b42 << ", b43 = " << b43 << std::endl;
}

int main() {
    // Parameters
    int N = 30;
    double tf = 10.0;
    double delta_v_max = 30.0;
    double delta_v_min = -30.0;
    double delta_s_max = 30.0;
    double delta_s_min = -30.0;
    double delta_m_max = 4000.0;
    double delta_m_min = -4000.0;
    double zmax = 0.0;
    double zmin = -100.0;
    double wmax = 5.0;
    double wmin = -5.0;
    double thetamax = 0.5;
    double thetamin = -0.5;
    double qmax = 5.0;
    double qmin = -5.0;
    double z0 = -20.0;
    double w0 = 0.1;
    double theta0 = 0.01;
    double q0 = 0.01;
    double delta_v0 = -20;
    double delta_s0 = -5.0;
    double delta_m0 = -300.00;
    double zf = -25.0;
    double thetaf = 0.0;
    double speed = 2.72;

    // Get state space matrices based on speed
    std::vector<double> A_flat, B_flat;
    getStateSpaceMatrices(speed, A_flat, B_flat);

    // Extract and print matrix elements
    extractAndPrintMatrixElements(A_flat, B_flat);

    return 0;
}
