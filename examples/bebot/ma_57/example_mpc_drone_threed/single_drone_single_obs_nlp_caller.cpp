#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

// ============================================================
// C API from libsingle_drone_single_obs_nlp.so
// ============================================================
extern "C" {
    struct PointSetProblem;

    PointSetProblem* create_point_set_problem(
        int N,
        double tf,

        double px_max,
        double px_min,
        double py_max,
        double py_min,
        double pz_max,
        double pz_min,
        double psi_max,
        double psi_min,

        double v_max,
        double w_max,
        double a_max,
        double aw_max,

        double px_cur,
        double py_cur,
        double pz_cur,
        double psi_cur,

        double vx_cur,
        double vy_cur,
        double vz_cur,
        double w_cur,

        double pxf,
        double pyf,
        double pzf,
        double psif,

        double sphere_x,
        double sphere_y,
        double sphere_z,
        double sphere_radius
    );

    int solve_point_set_problem(
        PointSetProblem* problem
    );

    int get_solution_size(
        PointSetProblem* problem
    );

    int get_solution(
        PointSetProblem* problem,
        double* solution,
        int n
    );

    double get_final_objective_value(
        PointSetProblem* problem
    );

    void destroy_point_set_problem(
        PointSetProblem* problem
    );
}

// ============================================================
// Bernstein helpers
// ============================================================
static std::vector<long double> binom_coeffs(int N) {
    std::vector<long double> C(N + 1, 0.0L);
    C[0] = 1.0L;

    for (int k = 1; k <= N; ++k) {
        C[k] =
            C[k - 1]
            * static_cast<long double>(N - k + 1)
            / static_cast<long double>(k);
    }

    return C;
}

static double bernstein_eval(
    const std::vector<double>& coeffs,
    double tau,
    double tau0,
    double tau1
) {
    const int degree =
        static_cast<int>(coeffs.size()) - 1;

    if (degree < 0) {
        return 0.0;
    }

    const double denom = tau1 - tau0;

    if (std::abs(denom) < 1e-12) {
        return coeffs.front();
    }

    double s =
        (tau - tau0) / denom;

    s = std::max(
        0.0,
        std::min(1.0, s)
    );

    const double one_minus_s =
        1.0 - s;

    const auto C =
        binom_coeffs(degree);

    long double value = 0.0L;

    for (int i = 0; i <= degree; ++i) {
        long double term =
            static_cast<long double>(coeffs[i])
            * C[i];

        term *= std::pow(
            static_cast<long double>(s),
            i
        );

        term *= std::pow(
            static_cast<long double>(one_minus_s),
            degree - i
        );

        value += term;
    }

    return static_cast<double>(value);
}

static bool extract_block(
    const std::vector<double>& solution,
    int block_index,
    int block_size,
    std::vector<double>& block
) {
    const int begin =
        block_index * block_size;

    const int end =
        begin + block_size;

    if (
        begin < 0
        || end > static_cast<int>(solution.size())
    ) {
        return false;
    }

    block.assign(
        solution.begin() + begin,
        solution.begin() + end
    );

    return true;
}

static double distance_3d(
    double x1,
    double y1,
    double z1,
    double x2,
    double y2,
    double z2
) {
    const double dx = x1 - x2;
    const double dy = y1 - y2;
    const double dz = z1 - z2;

    return std::sqrt(
        dx * dx
        + dy * dy
        + dz * dz
    );
}

// ============================================================
// Main
// ============================================================
int main() {
    // --------------------------------------------------------
    // OCP polynomial order / final transformed time
    // --------------------------------------------------------
    const int N = 4;
    const double tf = 0.9;

    const int L = N + 1;
    const int expected_solution_size = 12 * L;

    // --------------------------------------------------------
    // State / control bounds
    // --------------------------------------------------------
    const double px_min =
        -std::numeric_limits<double>::infinity();

    const double px_max =
        std::numeric_limits<double>::infinity();

    const double py_min =
        -std::numeric_limits<double>::infinity();

    const double py_max =
        std::numeric_limits<double>::infinity();

    const double pz_min = 0.0;
    const double pz_max = 5.0;

    const double psi_min =
        -3.14159265358979323846;

    const double psi_max =
        3.14159265358979323846;

    const double v_max = 1.0;
    const double w_max = 2.0;

    const double a_max = 1.5;
    const double aw_max = 8.0;

    // ========================================================
    // ONE INITIAL CONDITION
    // ========================================================
    const double px_cur  = -1.2;
    const double py_cur  =  1.8;
    const double pz_cur  =  0.6;
    const double psi_cur =  0.0;

    const double vx_cur = 0.0;
    const double vy_cur = 0.0;
    const double vz_cur = 0.0;
    const double w_cur  = 0.0;

    // ========================================================
    // ONE TARGET
    // ========================================================
    const double pxf  = 1.2;
    const double pyf  = 1.8;
    const double pzf  = 0.6;
    const double psif = 0.0;

    // ========================================================
    // ONE SPHERICAL OBSTACLE
    //
    // The obstacle lies directly between start and target.
    // The straight-line path would pass through its center.
    // ========================================================
    const double sphere_x = 0.0;
    const double sphere_y = 1.8;
    const double sphere_z = 0.6;
    const double sphere_radius = 0.35;

    // --------------------------------------------------------
    // Output settings
    //
    // Save the full optimized trajectory from tau = -1 to 0.9.
    // --------------------------------------------------------
    const int trajectory_samples = 1000;

    const double tau_start = -1.0;
    const double tau_end   =  0.9;

    const std::string output_csv =
        "single_drone_single_obs_nlp_trajectory.csv";

    // --------------------------------------------------------
    // Print test definition
    // --------------------------------------------------------
    std::cout
        << std::fixed
        << std::setprecision(10);

    std::cout
        << "============================================================\n"
        << "Single OCP solve: one start, one target, one sphere\n"
        << "============================================================\n"
        << "N  = " << N << "\n"
        << "tf = " << tf << "\n\n"

        << "Initial position:\n"
        << "  px = " << px_cur << "\n"
        << "  py = " << py_cur << "\n"
        << "  pz = " << pz_cur << "\n"
        << "  psi = " << psi_cur << "\n\n"

        << "Target:\n"
        << "  x = " << pxf << "\n"
        << "  y = " << pyf << "\n"
        << "  z = " << pzf << "\n"
        << "  psi = " << psif << "\n\n"

        << "Sphere obstacle:\n"
        << "  x = " << sphere_x << "\n"
        << "  y = " << sphere_y << "\n"
        << "  z = " << sphere_z << "\n"
        << "  radius = " << sphere_radius << "\n"
        << "============================================================\n\n";

    // ========================================================
    // CREATE ONE OCP
    // ========================================================
    PointSetProblem* problem =
        create_point_set_problem(
            N,
            tf,

            px_max,
            px_min,
            py_max,
            py_min,
            pz_max,
            pz_min,
            psi_max,
            psi_min,

            v_max,
            w_max,
            a_max,
            aw_max,

            px_cur,
            py_cur,
            pz_cur,
            psi_cur,

            vx_cur,
            vy_cur,
            vz_cur,
            w_cur,

            pxf,
            pyf,
            pzf,
            psif,

            sphere_x,
            sphere_y,
            sphere_z,
            sphere_radius
        );

    if (!problem) {
        std::cerr
            << "ERROR: create_point_set_problem() returned nullptr.\n";

        return 1;
    }

    // ========================================================
    // SOLVE ONCE
    // ========================================================
    const int solve_ok =
        solve_point_set_problem(problem);

    if (!solve_ok) {
        std::cerr
            << "ERROR: IPOPT OCP solve failed.\n";

        destroy_point_set_problem(problem);
        problem = nullptr;

        return 2;
    }

    // ========================================================
    // GET SOLUTION
    // ========================================================
    const int solution_size =
        get_solution_size(problem);

    if (solution_size != expected_solution_size) {
        std::cerr
            << "ERROR: unexpected solution size.\n"
            << "Expected: "
            << expected_solution_size
            << "\n"
            << "Received: "
            << solution_size
            << "\n";

        destroy_point_set_problem(problem);
        problem = nullptr;

        return 3;
    }

    std::vector<double> solution(
        solution_size,
        0.0
    );

    const int copied =
        get_solution(
            problem,
            solution.data(),
            solution_size
        );

    if (copied != solution_size) {
        std::cerr
            << "ERROR: could not copy full solution.\n";

        destroy_point_set_problem(problem);
        problem = nullptr;

        return 4;
    }

    const double objective =
        get_final_objective_value(problem);

    // ========================================================
    // EXTRACT THE 12 BERNSTEIN COEFFICIENT BLOCKS
    // ========================================================
    std::vector<double> px;
    std::vector<double> py;
    std::vector<double> pz;
    std::vector<double> psi;

    std::vector<double> vx;
    std::vector<double> vy;
    std::vector<double> vz;
    std::vector<double> w;

    std::vector<double> ax;
    std::vector<double> ay;
    std::vector<double> az;
    std::vector<double> aw;

    const bool blocks_ok =
        extract_block(solution,  0, L, px)
        && extract_block(solution,  1, L, py)
        && extract_block(solution,  2, L, pz)
        && extract_block(solution,  3, L, psi)

        && extract_block(solution,  4, L, vx)
        && extract_block(solution,  5, L, vy)
        && extract_block(solution,  6, L, vz)
        && extract_block(solution,  7, L, w)

        && extract_block(solution,  8, L, ax)
        && extract_block(solution,  9, L, ay)
        && extract_block(solution, 10, L, az)
        && extract_block(solution, 11, L, aw);

    if (!blocks_ok) {
        std::cerr
            << "ERROR: failed to extract one or more coefficient blocks.\n";

        destroy_point_set_problem(problem);
        problem = nullptr;

        return 5;
    }

    // --------------------------------------------------------
    // Print the position coefficient vectors.
    // --------------------------------------------------------
    auto print_coeffs =
        [](const std::string& name,
           const std::vector<double>& coeffs) {

            std::cout << name << " = [";

            for (size_t i = 0; i < coeffs.size(); ++i) {
                std::cout << coeffs[i];

                if (i + 1 < coeffs.size()) {
                    std::cout << ", ";
                }
            }

            std::cout << "]\n";
        };

    std::cout
        << "\nOptimal Bernstein position coefficients:\n";

    print_coeffs("px", px);
    print_coeffs("py", py);
    print_coeffs("pz", pz);

    std::cout
        << "\nObjective = "
        << objective
        << "\n";

    // ========================================================
    // SAVE THE FULL OPTIMIZED TRAJECTORY
    // ========================================================
    std::ofstream csv(output_csv);

    if (!csv.is_open()) {
        std::cerr
            << "ERROR: failed to open output CSV: "
            << output_csv
            << "\n";

        destroy_point_set_problem(problem);
        problem = nullptr;

        return 6;
    }

    csv
        << std::fixed
        << std::setprecision(10);

    csv
        << "iteration,"
        << "segment_sample,"
        << "tau,"
        << "x,"
        << "y,"
        << "z,"
        << "psi,"
        << "vx,"
        << "vy,"
        << "vz,"
        << "w,"
        << "ax,"
        << "ay,"
        << "az,"
        << "aw,"
        << "target_x,"
        << "target_y,"
        << "target_z,"
        << "target_psi,"
        << "sphere_x,"
        << "sphere_y,"
        << "sphere_z,"
        << "sphere_radius,"
        << "objective,"
        << "full_ocp_min_sphere_distance\n";

    double min_sphere_distance =
        std::numeric_limits<double>::infinity();

    // First evaluate all points so we know the minimum distance.
    for (int k = 0; k < trajectory_samples; ++k) {
        const double alpha =
            static_cast<double>(k)
            / static_cast<double>(
                trajectory_samples - 1
            );

        const double tau =
            tau_start
            + alpha * (tau_end - tau_start);

        const double x_eval =
            bernstein_eval(
                px,
                tau,
                -1.0,
                1.0
            );

        const double y_eval =
            bernstein_eval(
                py,
                tau,
                -1.0,
                1.0
            );

        const double z_eval =
            bernstein_eval(
                pz,
                tau,
                -1.0,
                1.0
            );

        const double distance =
            distance_3d(
                x_eval,
                y_eval,
                z_eval,
                sphere_x,
                sphere_y,
                sphere_z
            );

        min_sphere_distance =
            std::min(
                min_sphere_distance,
                distance
            );
    }

    // Now write the complete trajectory.
    for (int k = 0; k < trajectory_samples; ++k) {
        const double alpha =
            static_cast<double>(k)
            / static_cast<double>(
                trajectory_samples - 1
            );

        const double tau =
            tau_start
            + alpha * (tau_end - tau_start);

        const double x_eval =
            bernstein_eval(
                px,
                tau,
                -1.0,
                1.0
            );

        const double y_eval =
            bernstein_eval(
                py,
                tau,
                -1.0,
                1.0
            );

        const double z_eval =
            bernstein_eval(
                pz,
                tau,
                -1.0,
                1.0
            );

        const double psi_eval =
            bernstein_eval(
                psi,
                tau,
                -1.0,
                1.0
            );

        const double vx_eval =
            bernstein_eval(
                vx,
                tau,
                -1.0,
                1.0
            );

        const double vy_eval =
            bernstein_eval(
                vy,
                tau,
                -1.0,
                1.0
            );

        const double vz_eval =
            bernstein_eval(
                vz,
                tau,
                -1.0,
                1.0
            );

        const double w_eval =
            bernstein_eval(
                w,
                tau,
                -1.0,
                1.0
            );

        const double ax_eval =
            bernstein_eval(
                ax,
                tau,
                -1.0,
                1.0
            );

        const double ay_eval =
            bernstein_eval(
                ay,
                tau,
                -1.0,
                1.0
            );

        const double az_eval =
            bernstein_eval(
                az,
                tau,
                -1.0,
                1.0
            );

        const double aw_eval =
            bernstein_eval(
                aw,
                tau,
                -1.0,
                1.0
            );

        csv
            << 1 << ","
            << k << ","
            << tau << ","

            << x_eval << ","
            << y_eval << ","
            << z_eval << ","
            << psi_eval << ","

            << vx_eval << ","
            << vy_eval << ","
            << vz_eval << ","
            << w_eval << ","

            << ax_eval << ","
            << ay_eval << ","
            << az_eval << ","
            << aw_eval << ","

            << pxf << ","
            << pyf << ","
            << pzf << ","
            << psif << ","

            << sphere_x << ","
            << sphere_y << ","
            << sphere_z << ","
            << sphere_radius << ","

            << objective << ","
            << min_sphere_distance
            << "\n";
    }

    csv.close();

    // ========================================================
    // REPORT RESULTS
    // ========================================================
    const double clearance =
        min_sphere_distance
        - sphere_radius;

    const double final_x =
        bernstein_eval(
            px,
            tau_end,
            -1.0,
            1.0
        );

    const double final_y =
        bernstein_eval(
            py,
            tau_end,
            -1.0,
            1.0
        );

    const double final_z =
        bernstein_eval(
            pz,
            tau_end,
            -1.0,
            1.0
        );

    const double final_target_distance =
        distance_3d(
            final_x,
            final_y,
            final_z,
            pxf,
            pyf,
            pzf
        );

    std::cout
        << "\n============================================================\n"
        << "ONE-SOLVE RESULT\n"
        << "============================================================\n"
        << "Minimum distance from optimized trajectory to sphere center:\n"
        << "  "
        << min_sphere_distance
        << " m\n\n"

        << "Sphere radius:\n"
        << "  "
        << sphere_radius
        << " m\n\n"

        << "Minimum clearance:\n"
        << "  "
        << clearance
        << " m\n\n"

        << "Trajectory point at tau = "
        << tau_end
        << ":\n"
        << "  x = " << final_x << "\n"
        << "  y = " << final_y << "\n"
        << "  z = " << final_z << "\n\n"

        << "Distance from final optimized point to target:\n"
        << "  "
        << final_target_distance
        << " m\n\n"

        << "Saved full optimized trajectory to:\n"
        << "  "
        << output_csv
        << "\n"
        << "============================================================\n";

    if (
        min_sphere_distance
        + 5e-3
        < sphere_radius
    ) {
        std::cerr
            << "\nWARNING: dense numerical verification found "
            << "sphere penetration greater than 0.005 m.\n";
    }

    destroy_point_set_problem(problem);
    problem = nullptr;

    return 0;
}