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
        double zmax, double zmin,
        double wmax, double wmin,
        double thetamax, double thetamin,
        double qmax, double qmin,
        double ymax, double ymin,
        double psimax, double psimin,
        double vmax, double vmin,
        double rmax, double rmin,
        double z0, double w0, double theta0, double q0,
        double y0, double psi0, double v0, double r0,
        double delta_v0, double delta_s0, double delta_m0, double delta_h0,
        double zf, double thetaf, double yf, double psif,
        // then 64 A‐entries:
        double a11,double a12,double a13,double a14,double a15,double a16,double a17,double a18,
        double a21,double a22,double a23,double a24,double a25,double a26,double a27,double a28,
        double a31,double a32,double a33,double a34,double a35,double a36,double a37,double a38,
        double a41,double a42,double a43,double a44,double a45,double a46,double a47,double a48,
        double a51,double a52,double a53,double a54,double a55,double a56,double a57,double a58,
        double a61,double a62,double a63,double a64,double a65,double a66,double a67,double a68,
        double a71,double a72,double a73,double a74,double a75,double a76,double a77,double a78,
        double a81,double a82,double a83,double a84,double a85,double a86,double a87,double a88,
        // then 32 B‐entries:
        double b11,double b12,double b13,double b14,
        double b21,double b22,double b23,double b24,
        double b31,double b32,double b33,double b34,
        double b41,double b42,double b43,double b44,
        double b51,double b52,double b53,double b54,
        double b61,double b62,double b63,double b64,
        double b71,double b72,double b73,double b74,
        double b81,double b82,double b83,double b84,
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
    double tf = 100.0;
    double t0 =  0.0, tend = tf;

    double delta_v_max =  30.0, delta_v_min = -30.0;
    double delta_s_max =  30.0, delta_s_min = -30.0;
    double delta_m_max =4000.0, delta_m_min=-4000.0;
    double delta_h_max =  30.0, delta_h_min = -30.0;

    double zmax=   0.0, zmin= -100.0;
    double wmax=   5.0, wmin=   -5.0;
    double thetamax=0.5, thetamin=-0.5;
    double qmax=   5.0, qmin=   -5.0;
    double ymax= 100.0, ymin=-100.0;
    double psimax=1.57, psimin=-1.57;
    double vmax=   5.0, vmin=   -5.0;
    double rmax=   5.0, rmin=   -5.0;

    double z0 = -20.0, w0 = -0.01, theta0 =  0.0,   q0 = 0.0;
    double y0 =   0.5, psi0 =  0.1,  v0 =  0.0,    r0 = 0.0;
    double dv0=   0.0, ds0 =  0.0,  dm0 =  0.0,  dh0 = 0.0;
    double zf = -25.0, thetaf = 0.0, yf = 0.0,   psif = 0.0;

    // A‐matrix (row‐major):
    double A[8][8] = {
        {-7.53007450746585e-15, 4.99999999991616,    0.999999999999994,  -3.65319419902939e-13, -9.99602909231959e-16, 4.31737682959773e-14, 2.66559016936731e-14, -2.68064244536696e-12},
        {-4.02549191160255e-17, 9.80614054063698e-15,2.28742631772763e-15, 0.999999999999981,   -1.51352208268968e-17, 3.33117779407292e-16, 4.95404385666044e-16,  2.07007768210703e-14},
        {0.000840272436887615, -3.04088623261381e-07,-0.101668320745939, -2.70007766110086,   -2.27750116762283e-17, -3.37167414420878e-16,-0.0656192217604660, -0.233946610893178},
        {4.37946882422992e-06,  -0.0149473707777599, -0.00572338624834895,-0.244530610026609,  -2.97838974510285e-18, 6.89414872980893e-17, -6.56780683930467e-06,-0.000331476363457423},
        {-5.01843558958588e-16,-9.98908090991322e-15,-1.99982289648071e-15,-6.15750233160021e-14,-1.24632263746486e-14, -4.99999999991671, 0.999999999999934,   3.70250104011795e-13},
        {2.31463391150403e-17, 2.00837561476649e-15, 6.09344582110067e-16, -3.07201372173251e-15, -5.52605011017448e-17, 2.52149400093135e-16, 2.29367009667294e-15,  1.00000000000006},
        {-1.71125752355647e-07, 8.11075295919674e-10, -0.00112476442426664,  0.0358936866747503, -5.36322241557343e-16, 4.53241816028819e-15, -0.0229753866801526,  6.37645253694799},
        {-2.50799373371796e-08, 9.15667953463273e-10, -0.000166140206403180, 0.00530192100279736, 1.41825673657912e-18, 7.23410560380495e-17, -0.00338810792080797, -0.158037916785118}
    };

    // B‐matrix (row‐major):
    double B[8][4] = {
        {-7.31017259460342e-17,-2.30059386320077e-20, 2.41878268201211e-18, -9.85959255647513e-17},
        {3.54645421442446e-20,-3.24012983442396e-22,  3.25883927778635e-19,  1.01964593573660e-19},
        {-0.00261112738112122,-2.21991902295883e-06, 0.00164151244056022, -3.73673141014142e-07},
        {0.000330132352958294,7.45681088024612e-23,  7.27061635497237e-05,  5.80930190419673e-07},
        {-4.50724817990609e-18,3.31749779096055e-21, -1.85246922534076e-18,  2.32598339321264e-17},
        {3.92236857955907e-20, -6.38268997174661e-23, 2.03133234178711e-20,  1.38993811781273e-19},
        {-9.81552974016019e-05,-1.95123434146839e-22,-4.70209873220215e-17,  0.00200500289959419},
        {-1.44986285170095e-05,-2.89891272257610e-23, 1.30016091925287e-16,  0.000295671463475263}
    };

    // flatten and pass everything
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
        ymax, ymin,
        psimax, psimin,
        vmax, vmin,
        rmax, rmin,
        z0, w0, theta0, q0,
        y0, psi0, v0, r0,
        dv0, ds0, dm0, dh0,
        zf, thetaf, yf, psif,

        // A entries
        A[0][0],A[0][1],A[0][2],A[0][3],A[0][4],A[0][5],A[0][6],A[0][7],
        A[1][0],A[1][1],A[1][2],A[1][3],A[1][4],A[1][5],A[1][6],A[1][7],
        A[2][0],A[2][1],A[2][2],A[2][3],A[2][4],A[2][5],A[2][6],A[2][7],
        A[3][0],A[3][1],A[3][2],A[3][3],A[3][4],A[3][5],A[3][6],A[3][7],
        A[4][0],A[4][1],A[4][2],A[4][3],A[4][4],A[4][5],A[4][6],A[4][7],
        A[5][0],A[5][1],A[5][2],A[5][3],A[5][4],A[5][5],A[5][6],A[5][7],
        A[6][0],A[6][1],A[6][2],A[6][3],A[6][4],A[6][5],A[6][6],A[6][7],
        A[7][0],A[7][1],A[7][2],A[7][3],A[7][4],A[7][5],A[7][6],A[7][7],

        // B entries
        B[0][0],B[0][1],B[0][2],B[0][3],
        B[1][0],B[1][1],B[1][2],B[1][3],
        B[2][0],B[2][1],B[2][2],B[2][3],
        B[3][0],B[3][1],B[3][2],B[3][3],
        B[4][0],B[4][1],B[4][2],B[4][3],
        B[5][0],B[5][1],B[5][2],B[5][3],
        B[6][0],B[6][1],B[6][2],B[6][3],
        B[7][0],B[7][1],B[7][2],B[7][3],

        t0, tend
    );

    // solve
    solve_point_set_problem(prob);

    // print objective
    std::cout << "Final objective: "
              << get_final_objective_value(prob) << "\n";

    // get / print solution vector
    int n_vars = 12*(N+1);
    std::vector<double> sol(n_vars);
    //get_solution(prob, sol.data(), n_vars);
    //std::cout << "Solution (first 12*(N+1) vars):\n";
    //for(int i=0; i<n_vars; ++i){
    //  std::cout << sol[i] << ( (i%10==9) ? "\n" : " " );
    //}

    //destroy_point_set_problem(prob);
    return 0;
}
