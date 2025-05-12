#include <iostream>
#include <vector>
#include "state_space_matrices.h"

// Declare the function from the shared library
extern "C" void create_point_set_problem(int N, double tf, double delta_v_max, double delta_v_min, double delta_s_max, double delta_s_min,
                                        double delta_m_max, double delta_m_min, double zmax, double zmin, double wmax, double wmin,
                                        double thetamax, double thetamin, double qmax, double qmin, double z0, double w0, double theta0,
                                        double q0, double delta_v0, double delta_s0, double delta_m0, double zf, double thetaf,
                                        double a11, double a12, double a13, double a14, double a21, double a22, double a23, double a24,
                                        double a31, double a32, double a33, double a34, double a41, double a42, double a43, double a44,
                                        double b11, double b12, double b13, double b21, double b22, double b23, double b31, double b32,
                                        double b33, double b41, double b42, double b43, double t0, double tend, double psi0);

void printMatrixElements(const std::vector<double>& A, const std::vector<double>& B) {
    // Extracting matrix elements
    double a11 = A[0], a12 = A[1], a13 = A[2], a14 = A[3];
    double a21 = A[4], a22 = A[5], a23 = A[6], a24 = A[7];
    double a31 = A[8], a32 = A[9], a33 = A[10], a34 = A[11];
    double a41 = A[12], a42 = A[13], a43 = A[14], a44 = A[15];

    double b11 = B[0], b12 = B[1], b13 = B[2];
    double b21 = B[3], b22 = B[4], b23 = B[5];
    double b31 = B[6], b32 = B[7], b33 = B[8];
    double b41 = B[9], b42 = B[10], b43 = B[11];

    // Print A elements
    std::cout << "Matrix A elements:" << std::endl;
    std::cout << "a11 = " << a11 << ", a12 = " << a12 << ", a13 = " << a13 << ", a14 = " << a14 << std::endl;
    std::cout << "a21 = " << a21 << ", a22 = " << a22 << ", a23 = " << a23 << ", a24 = " << a24 << std::endl;
    std::cout << "a31 = " << a31 << ", a32 = " << a32 << ", a33 = " << a33 << ", a34 = " << a34 << std::endl;
    std::cout << "a41 = " << a41 << ", a42 = " << a42 << ", a43 = " << a43 << ", a44 = " << a44 << std::endl;

    // Print B elements
    std::cout << "Matrix B elements:" << std::endl;
    std::cout << "b11 = " << b11 << ", b12 = " << b12 << ", b13 = " << b13 << std::endl;
    std::cout << "b21 = " << b21 << ", b22 = " << b22 << ", b23 = " << b23 << std::endl;
    std::cout << "b31 = " << b31 << ", b32 = " << b32 << ", b33 = " << b33 << std::endl;
    std::cout << "b41 = " << b41 << ", b42 = " << b42 << ", b43 = " << b43 << std::endl;
}

int main() {
    int N = 20;
    double tf = 25.0;
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
    double t0 = 0;
    double tend = 10.0;
    double psi0 = 0;

    std::vector<double> A_flat, B_flat;
    getStateSpaceMatrices(speed, A_flat, B_flat);

    double a11 = A_flat[0], a12 = A_flat[1], a13 = A_flat[2], a14 = A_flat[3];
    double a21 = A_flat[4], a22 = A_flat[5], a23 = A_flat[6], a24 = A_flat[7];
    double a31 = A_flat[8], a32 = A_flat[9], a33 = A_flat[10], a34 = A_flat[11];
    double a41 = A_flat[12], a42 = A_flat[13], a43 = A_flat[14], a44 = A_flat[15];

    double b11 = B_flat[0], b12 = B_flat[1], b13 = B_flat[2];
    double b21 = B_flat[3], b22 = B_flat[4], b23 = B_flat[5];
    double b31 = B_flat[6], b32 = B_flat[7], b33 = B_flat[8];
    double b41 = B_flat[9], b42 = B_flat[10], b43 = B_flat[11];

    create_point_set_problem(N, tf, delta_v_max, delta_v_min, delta_s_max, delta_s_min, delta_m_max, delta_m_min, zmax, zmin, wmax, wmin,
                            thetamax, thetamin, qmax, qmin, z0, w0, theta0, q0, delta_v0, delta_s0, delta_m0, zf, thetaf,
                            a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34, a41, a42, a43, a44,
                            b11, b12, b13, b21, b22, b23, b31, b32, b33, b41, b42, b43, t0, tend, psi0);

    return 0;
}
