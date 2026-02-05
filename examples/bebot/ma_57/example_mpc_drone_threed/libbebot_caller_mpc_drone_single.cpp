#include <iostream>
#include <vector>

// declare the C‐API from your shared library
extern "C" {
    struct PointSetProblem;

    PointSetProblem* create_point_set_problem(
        int N, double tf,
        double delta_v_max, double delta_v_min,
        double delta_s_max, double delta_s_min,
        double delta_m_max, double delta_m_min,
        double delta_h_max, double delta_h_min,
        // double delta_n_max, double delta_n_min,
        double zmax, double zmin,
        double wmax, double wmin,
        double thetamax, double thetamin,
        double qmax, double qmin,
        // double umax, double umin,
        double psimax, double psimin,
        double rmax, double rmin,
        // double xmax, double xmin, double ymax, double ymin,
        double z0, double w0, double theta0, double q0,
        // double u0,
        double psi0, double r0,
        // double x0, double y0,
        double delta_v0, double delta_s0, double delta_m0, double delta_h0,
        // double delta_n0,
        double zf, double thetaf,
        // double xf, double yf,
        double psif,

        // A (4x4)
        double a11,double a12,double a13,double a14,
        double a21,double a22,double a23,double a24,
        double a31,double a32,double a33,double a34,
        double a41,double a42,double a43,double a44,

        // B (4x3)
        double b11,double b12,double b13,
        double b21,double b22,double b23,
        double b31,double b32,double b33,
        double b41,double b42,double b43,

        // C (2x2)
        double c11,double c12,
        double c21,double c22,

        // D (2x1)
        double d11,
        double d21,

        // finally t0, tend,
        double t0, double tend,

        // ---- NEW: previous controls (fixed at first knot) ----
        double dv_prev, double ds_prev, double dm_prev, double dh_prev,

        // ---- NEW: control time-derivative bounds ----
        double dv_dot_max, double dm_dot_max, double ds_dot_max, double dh_dot_max
    );

    void solve_point_set_problem(PointSetProblem* problem);
    void get_solution(PointSetProblem* problem, double* solution, int n);
    double get_final_objective_value(PointSetProblem* problem);
    void destroy_point_set_problem(PointSetProblem* problem);
}

int main() {
    // problem dimensions & bounds
    int    N  = 5;
    double tf = 10.0;
    double t0 = 0.0, tend = 10.0;

    double delta_v_max =  30.0, delta_v_min = -30.0;
    double delta_s_max =  30.0, delta_s_min = -30.0;
    double delta_m_max = 4000.0, delta_m_min = -4000.0;
    double delta_h_max =  30.0, delta_h_min = -30.0;

    double zmax=   0.0, zmin= -100.0;
    double wmax=  10.0, wmin=  -10.0;
    double thetamax=0.5, thetamin=-0.5;
    double qmax=  10.0, qmin=  -10.0;

    double psimax= 3.14, psimin=-3.14;
    double rmax=   5.0,  rmin=  -5.0;

    double z0 = -19.97, w0 = 0.02, theta0 = 0.00, q0 = 0.0;
    double psi0 = 0.0, r0 = 0.0;

    // These "delta_*0" are constructor inputs, but your bounds actually pin knot-0 to *_prev (below).
    double dv0 = 0.0, ds0 = 0.0, dm0 = 0.0, dh0 = 0.0;

    double zf = -23.0, thetaf = 0.0;
    double psif = 0.0;

    // NEW: previous controls (used in get_bounds_info: x_lower[...]=x_upper[...]=*_prev)
    double dv_prev = 0.13;
    double ds_prev = 4.02;
    double dm_prev = -2.86;
    double dh_prev = 0.42;

    // NEW: derivative bounds (used as bounds on g blocks 6..9)
    double dv_dot_max = 0.6;
    double dm_dot_max = 8.5;
    double ds_dot_max = 0.6;
    double dh_dot_max = 0.6;

    // A‐matrix (row‐major):
    double A[4][4] = {
        {-9.61335633071093e-13,  4.99999999991687,   1.00000000005217,  1.12249207912694e-10},
        { 1.45212601715332e-12, -5.13981176158834e-14,-1.17800813245335e-10, 1.00000000101190},
        { 0.000840272043916469, -3.04085570568224e-07,-0.101668288331192, -2.70007795164826},
        { 4.37929332288747e-06, -0.0149473707639360, -0.00572337288021148,-0.244530705224609}
    };

    // B‐matrix (row‐major):
    double B[4][3] = {
        {-1.69674752574204e-12, -1.10484134201392e-17,  2.75706811981441e-12},
        {-7.22760359972355e-13, -6.12418309548146e-18, -4.20492416927654e-12},
        {-0.00261112713870843,  -2.21991902099372e-06,  0.00164151357821695},
        { 0.000330132364452263,  2.40305021375118e-16,  7.27066721427474e-05}
    };

    // C‐matrix (row‐major):
    double C[2][2] = {
        {-3.44776290850390e-16,  1.00000000000000},
        { 9.60275984420308e-17, -0.158037916785125}
    };

    // D‐matrix (row‐major):
    double D[2][1] = {
        {1.35525271560688e-19},
        {0.000295671463475263}
    };

    // create problem
    PointSetProblem* prob = create_point_set_problem(
        N, tf,
        delta_v_max, delta_v_min,
        delta_s_max, delta_s_min,
        delta_m_max, delta_m_min,
        delta_h_max, delta_h_min,

        zmax, zmin,
        wmax, wmin,
        thetamax, thetamin,
        qmax, qmin,

        psimax, psimin,
        rmax, rmin,

        z0, w0, theta0, q0,
        psi0, r0,

        dv0, ds0, dm0, dh0,

        zf, thetaf,
        psif,

        // A entries
        A[0][0],A[0][1],A[0][2],A[0][3],
        A[1][0],A[1][1],A[1][2],A[1][3],
        A[2][0],A[2][1],A[2][2],A[2][3],
        A[3][0],A[3][1],A[3][2],A[3][3],

        // B entries
        B[0][0],B[0][1],B[0][2],
        B[1][0],B[1][1],B[1][2],
        B[2][0],B[2][1],B[2][2],
        B[3][0],B[3][1],B[3][2],

        // C entries
        C[0][0],C[0][1],
        C[1][0],C[1][1],

        // D entries
        D[0][0],
        D[1][0],

        // time window
        t0, tend,

        // NEW: previous controls
        dv_prev, ds_prev, dm_prev, dh_prev,

        // NEW: derivative bounds
        dv_dot_max, dm_dot_max, ds_dot_max, dh_dot_max
    );

    // solve
    solve_point_set_problem(prob);

    // print objective
    std::cout << "Final objective: "
              << get_final_objective_value(prob) << "\n";

    // get / print solution vector
    // int n_vars = 10*(N+1);
    // std::vector<double> sol(n_vars);
    //get_solution(prob, sol.data(), n_vars);
    //std::cout << "Solution (first 12*(N+1) vars):\n";
    //for(int i=0; i<n_vars; ++i){
    //  std::cout << sol[i] << ( (i%10==9) ? "\n" : " " );
    //}

    //destroy_point_set_problem(prob);
    return 0;
}
