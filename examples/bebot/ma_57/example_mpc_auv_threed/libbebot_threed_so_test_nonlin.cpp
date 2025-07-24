// main.cpp
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
        double delta_n_max, double delta_n_min,
        double zmax, double zmin,
        double wmax, double wmin,
        double thetamax, double thetamin,
        double qmax, double qmin,
        double umax, double umin,
        double psimax, double psimin,
        double rmax, double rmin,
        double xmax, double xmin, double ymax, double ymin,
        double z0, double w0, double theta0, double q0,
        double u0, double psi0, double r0,
        double x0, double y0, 
        double delta_v0, double delta_s0, double delta_m0, double delta_h0, double delta_n0,
        double zf, double thetaf, double xf, double yf, double psif,
        // then 64 A‐entries:
        double a11,double a12,double a13,double a14,
        double a21,double a22,double a23,double a24,
        double a31,double a32,double a33,double a34,
        double a41,double a42,double a43,double a44,
        // then 32 B‐entries:
        double b11,double b12,double b13,
        double b21,double b22,double b23,
        double b31,double b32,double b33,
        double b41,double b42,double b43,
        // C
        double c11,double c12,double c13,
        double c21,double c22,double c23,
        double c31,double c32,double c33,
        // d
        double d11,double d12,
        double d21,double d22,
        double d31,double d32,
        
        // finally t0, tend
        double t0, double tend
    );
    void solve_point_set_problem(PointSetProblem* problem);
    void get_solution(PointSetProblem* problem, double* solution, int n);
    double get_final_objective_value(PointSetProblem* problem);
    void destroy_point_set_problem(PointSetProblem* problem);
}

int main() {
    // problem dimensions & bounds
    int    N  = 10;
    double tf = 200.0;
    double t0 =  0.0, tend = tf;

    double delta_v_max =  30.0, delta_v_min = -30.0;
    double delta_s_max =  30.0, delta_s_min = -30.0;
    double delta_m_max =4000.0, delta_m_min=-4000.0;
    double delta_h_max =  30.0, delta_h_min = -30.0;
    double delta_n_max =  -1.1, delta_n_min = -1.3;

    double zmax=   0.0, zmin= -100.0;
    double wmax=   5.0, wmin=   -5.0;
    double thetamax=0.5, thetamin=-0.5;
    double qmax=   5.0, qmin=   -5.0;

    double umax=   -5.1, umin=   -5.5;
    double psimax=3.14, psimin=-3.14;
    double rmax=   1.0, rmin=   -1.0;

    double xmax= 5000, xmin=-5000.0;
    double ymax= 100.0, ymin=-100.0;

    double z0 = -25.00, w0 = 0.00, theta0 =  0.0,   q0 = 0.0;
    
    double u0 =  5.115, psi0 =  0.0,  r0 = 0.0;
    double x0 =   0.0, y0 =   0.0;
    
    double dv0=   0.0, ds0 =  0.0,  dm0 =  0.0,  dh0 = 0.0, dn0 = -1.25;
    double zf = -20.0, thetaf = 0.0, xf = -1000.0, yf = 00.0,   psif = 0.0;

    // A‐matrix (row‐major):
    double A[4][4] = {
        {-9.61335633071093e-13,	4.99999999991687,   1.00000000005217,	1.12249207912694e-10},
        {1.45212601715332e-12,    -5.13981176158834e-14,	-1.17800813245335e-10,	1.00000000101190},
        {0.000840272043916469,    -3.04085570568224e-07,	-0.101668288331192,	-2.70007795164826},
        {4.37929332288747e-06,    -0.0149473707639360,	-0.00572337288021148,	-0.244530705224609}
    };

    // B‐matrix (row‐major):
    double B[4][3] = {
        {-1.69674752574204e-12,	-1.10484134201392e-17,	2.75706811981441e-12},
        {-7.22760359972355e-13,	-6.12418309548146e-18,	-4.20492416927654e-12},
        {-0.00261112713870843,     -2.21991902099372e-06,	0.00164151357821695},
        {0.000330132364452263,	2.40305021375118e-16,	7.27066721427474e-05}
    };

    // C‐matrix (row‐major):
    double C[3][3] = {
        {-0.00940468274749529,	1.90475123646174e-15,	6.15646652556745e-05},
        {1.91915088359914e-15,	-7.30423863335897e-16,	1.00000000000003},
        {-3.32911359041052e-07,	1.00161193136768e-16,	-0.158037916785119}
    };

    // D‐matrix (row‐major):
    double D[3][2] = {
        {6.60883567194786e-07,	0.0386000677170305},
        {2.18137750750520e-19,	1.87097290348107e-18},
        {0.000295671463475263,	-3.19320363176623e-19}
    };


    // flatten and pass everything
    PointSetProblem* prob = create_point_set_problem(
        N, tf,
        delta_v_max, delta_v_min,
        delta_s_max, delta_s_min,
        delta_m_max, delta_m_min,
        delta_h_max, delta_h_min,
        delta_n_max, delta_n_min,
        zmax, zmin,
        wmax, wmin,
        thetamax, thetamin,
        qmax, qmin,
        umax, umin,
        psimax, psimin,
        rmax, rmin,
        xmax, xmin,
        ymax, ymin,
        z0, w0, theta0, q0,
        u0, psi0, r0,
        x0, y0, 
        dv0, ds0, dm0, dh0, dn0,
        zf, thetaf, xf, yf, psif,

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
        C[0][0],C[0][1],C[0][2],
        C[1][0],C[1][1],C[1][2],
        C[2][0],C[2][1],C[2][2],

        // D entries
        D[0][0],D[0][1],
        D[1][0],D[1][1],
        D[2][0],D[2][1],
        
        
        t0, tend
    );

    // solve
    solve_point_set_problem(prob);

    // print objective
    std::cout << "Final objective: "
              << get_final_objective_value(prob) << "\n";

    // get / print solution vector
    int n_vars = 14*(N+1) + 1;
    std::vector<double> sol(n_vars);
    //get_solution(prob, sol.data(), n_vars);
    //std::cout << "Solution (first 12*(N+1) vars):\n";
    //for(int i=0; i<n_vars; ++i){
    //  std::cout << sol[i] << ( (i%10==9) ? "\n" : " " );
    //}

    //destroy_point_set_problem(prob);
    return 0;
}
