#include <iostream>
#include <vector>
#include <iomanip>
#include <limits>

#include <fstream>

#include <cmath>

// C-API from your shared library
extern "C" {
    struct PointSetProblem;

    PointSetProblem* create_point_set_problem(
        int N, double tf,
        double px_max, double px_min,
        double py_max, double py_min,
        double pz_max, double pz_min,
        double psi_max, double psi_min,
        double v_max,  double w_max,
        double a_max,  double aw_max,
        double px_cur, double py_cur, double pz_cur, double psi_cur,
        double vx_cur, double vy_cur, double vz_cur, double w_cur,
        double pxf,    double pyf,    double pzf,    double psif
    );

    void   solve_point_set_problem(PointSetProblem* problem);
    void   get_solution(PointSetProblem* problem, double* solution, int n);
    double get_final_objective_value(PointSetProblem* problem);
    void   destroy_point_set_problem(PointSetProblem* problem);
}

int main() {
    // -----------------------------
    // Problem size
    // -----------------------------
    const int    N  = 4;
    const double tf = 0.9;

    // -----------------------------
    // Bounds (example values — set to what YOU want)
    // -----------------------------
    const double px_min = -std::numeric_limits<double>::infinity(),  px_max =  std::numeric_limits<double>::infinity();
    const double py_min = -std::numeric_limits<double>::infinity(),  py_max =  std::numeric_limits<double>::infinity();
    const double pz_min = 0,  pz_max =  5.0;

    const double psi_min = -3.141592653589793, psi_max = 3.141592653589793;

    const double v_max  = 1.0;   // bounds vx,vy,vz in [-v_max, +v_max]
    const double w_max  = 2.0;   // bounds w in [-w_max, +w_max]

    const double a_max  = 1.5;   // bounds ax,ay,az in [-a_max, +a_max]
    const double aw_max = 8.0;   // bounds aw in [-aw_max, +aw_max]

    // -----------------------------
    // Current state at knot 0 (these are pinned in your get_bounds_info)
    // -----------------------------
    const double px_cur  = 0.2023929884; // 		
    const double py_cur  = 1.133773398;
    const double pz_cur  = 0.6;
    const double psi_cur = 0.0;

    const double vx_cur  = 0.0;
    const double vy_cur  = 0.0;
    const double vz_cur  = 0.0;
    const double w_cur   = 0.0;

    // -----------------------------
    // Target (used in objective)
    // -----------------------------
    const double pxf  = 1.0;
    const double pyf  = 2.5;
    const double pzf  = 0.6;
    const double psif = 0.0;

    // -----------------------------
    // Create + solve
    // -----------------------------
    PointSetProblem* prob = create_point_set_problem(
        N, tf,
        px_max, px_min,
        py_max, py_min,
        pz_max, pz_min,
        psi_max, psi_min,
        v_max, w_max,
        a_max, aw_max,
        px_cur, py_cur, pz_cur, psi_cur,
        vx_cur, vy_cur, vz_cur, w_cur,
        pxf, pyf, pzf, psif
    );

    if (!prob) {
        std::cerr << "create_point_set_problem() returned nullptr\n";
        return 1;
    }

    solve_point_set_problem(prob);

    const double J = get_final_objective_value(prob);
    std::cout << std::setprecision(16)
              << "Final objective: " << J << "\n";

    // -----------------------------
    // Fetch solution vector (12*(N+1))
    // Layout in your TNLP: [px,py,pz,psi,vx,vy,vz,w,ax,ay,az,aw] each length (N+1)
    // -----------------------------
    // const int n_vars = 12 * (N + 1);
    // std::vector<double> sol(n_vars, 0.0);

    // get_solution(prob, sol.data(), n_vars);

    // std::cout << "Solution control points (first block px):\n";
    // for (int i = 0; i < (N + 1); ++i) {
    //     std::cout << "px[" << i << "] = " << sol[i] << "\n";
    // }

    // destroy_point_set_problem(prob);
    return 0;
}


// g++ -O2 -std=c++17 -o single_drone_offline_test libbebot_caller_mpc_drone_single.cpp   -L. -lbebot_mpc_drone_single -Wl,-rpath,'$ORIGIN'
