// #include <iostream>
// #include <vector>
// #include <iomanip>
// #include <limits>
// #include <fstream>
// #include <cmath>
// #include <algorithm>
// #include <string>

// // C-API from your shared library
// extern "C" {
//     struct PointSetProblem;

//     PointSetProblem* create_point_set_problem(
//         int N, double tf,
//         double px_max, double px_min,
//         double py_max, double py_min,
//         double pz_max, double pz_min,
//         double psi_max, double psi_min,
//         double v_max,  double w_max,
//         double a_max,  double aw_max,
//         double px_cur, double py_cur, double pz_cur, double psi_cur,
//         double vx_cur, double vy_cur, double vz_cur, double w_cur,
//         double pxf,    double pyf,    double pzf,    double psif,

//         // Cylinder obstacle parameters
//         double cyl_x,
//         double cyl_y,
//         double cyl_radius,

//         // Sphere obstacle parameters
//         double sphere_x,
//         double sphere_y,
//         double sphere_z,
//         double sphere_radius
//     );

//     void   solve_point_set_problem(PointSetProblem* problem);
//     void   get_solution(PointSetProblem* problem, double* solution, int n);
//     double get_final_objective_value(PointSetProblem* problem);
//     void   destroy_point_set_problem(PointSetProblem* problem);
// }

// // ------------------------------------------------------------
// // Scenario:
// // target + one cylinder slot + one sphere slot
// // ------------------------------------------------------------
// struct MpcScenario {
//     double pxf;
//     double pyf;
//     double pzf;
//     double psif;

//     double cyl_x;
//     double cyl_y;
//     double cyl_radius;

//     double sphere_x;
//     double sphere_y;
//     double sphere_z;
//     double sphere_radius;
// };

// // ------------------------------------------------------------
// // Safe Bernstein evaluator
// // ------------------------------------------------------------
// static std::vector<long double> binom_coeffs(int N) {
//     std::vector<long double> C(N + 1, 0.0L);
//     C[0] = 1.0L;

//     for (int k = 1; k <= N; ++k) {
//         C[k] = C[k - 1]
//              * static_cast<long double>(N - k + 1)
//              / static_cast<long double>(k);
//     }

//     return C;
// }

// static double bernstein_eval(
//     const std::vector<double>& coeffs,
//     double t,
//     double t0,
//     double tf
// ) {
//     const int degree = static_cast<int>(coeffs.size()) - 1;

//     if (degree < 0) {
//         return 0.0;
//     }

//     const double denom = tf - t0;

//     if (std::abs(denom) < 1e-12) {
//         return coeffs.front();
//     }

//     double s = (t - t0) / denom;

//     if (s < 0.0) {
//         s = 0.0;
//     }

//     if (s > 1.0) {
//         s = 1.0;
//     }

//     const double one_minus_s = 1.0 - s;
//     const auto C = binom_coeffs(degree);

//     long double value = 0.0L;

//     for (int i = 0; i <= degree; ++i) {
//         long double term = static_cast<long double>(coeffs[i]) * C[i];

//         term *= ::powl(
//             static_cast<long double>(s),
//             static_cast<long double>(i)
//         );

//         term *= ::powl(
//             static_cast<long double>(one_minus_s),
//             static_cast<long double>(degree - i)
//         );

//         value += term;
//     }

//     return static_cast<double>(value);
// }

// // ------------------------------------------------------------
// // Extract one block from solution.
// // Layout:
// //   block 0: px
// //   block 1: py
// //   block 2: pz
// //   block 3: psi
// //   block 4: vx
// //   block 5: vy
// //   block 6: vz
// //   block 7: w
// //   block 8: ax
// //   block 9: ay
// //   block 10: az
// //   block 11: aw
// // ------------------------------------------------------------
// static bool extract_block(
//     const std::vector<double>& sol,
//     int block_idx,
//     int Np1,
//     std::vector<double>& out
// ) {
//     const int start = block_idx * Np1;
//     const int end   = (block_idx + 1) * Np1;

//     if (start < 0 || end > static_cast<int>(sol.size()) || end <= start) {
//         return false;
//     }

//     out.assign(sol.begin() + start, sol.begin() + end);
//     return true;
// }

// static void print_vector(
//     const std::string& name,
//     const std::vector<double>& v
// ) {
//     std::cout << name << " = [";

//     for (size_t i = 0; i < v.size(); ++i) {
//         std::cout << std::setprecision(10) << v[i];

//         if (i + 1 < v.size()) {
//             std::cout << ", ";
//         }
//     }

//     std::cout << "]\n";
// }

// int main() {
//     // ------------------------------------------------------------
//     // Receding-horizon loop settings
//     //
//     // Scenario 1:
//     //   obstacle set 1 + target 1
//     //
//     // Scenario 2:
//     //   obstacle set 2 + target 2
//     //
//     // Scenario 3:
//     //   no obstacles + same target as scenario 2
//     // ------------------------------------------------------------
//     const int scenario1_iterations = 3;
//     const int scenario2_iterations = 4;
//     const int scenario3_iterations = 15;

//     const int num_mpc_iterations =
//         scenario1_iterations
//       + scenario2_iterations
//       + scenario3_iterations;

//     const double Ts = 0.05;
//     const double sample_index = 24.0;
//     const double t_sample = (sample_index - 1.0) * Ts;
//     const double tau_sample = (t_sample - 1.0) / (t_sample + 1.0);

//     const bool update_velocity_state_from_solution = true;

//     // Optional:
//     // If true, velocity is reset to zero exactly when scenario changes.
//     // Usually false is more dynamically consistent.
//     const bool reset_velocity_when_scenario_changes = false;

//     // ------------------------------------------------------------
//     // Problem size
//     // ------------------------------------------------------------
//     const int    N  = 4;
//     const double tf = 0.9;

//     const int Np1 = N + 1;
//     const int solution_size = 12 * Np1;

//     // ------------------------------------------------------------
//     // Bounds
//     // ------------------------------------------------------------
//     const double px_min = -std::numeric_limits<double>::infinity();
//     const double px_max =  std::numeric_limits<double>::infinity();

//     const double py_min = -std::numeric_limits<double>::infinity();
//     const double py_max =  std::numeric_limits<double>::infinity();

//     const double pz_min = 0.0;
//     const double pz_max = 5.0;

//     const double psi_min = -3.141592653589793;
//     const double psi_max =  3.141592653589793;

//     const double v_max  = 1.0;
//     const double w_max  = 2.0;

//     const double a_max  = 1.5;
//     const double aw_max = 8.0;

//     // ------------------------------------------------------------
//     // Initial state for first OCP
//     // ------------------------------------------------------------
//     double px_cur  = -1.2;
//     double py_cur  =  1.8;
//     double pz_cur  =  0.6;
//     double psi_cur =  0.0;

//     double vx_cur = 0.0;
//     double vy_cur = 0.0;
//     double vz_cur = 0.0;
//     double w_cur  = 0.0;

//     // ------------------------------------------------------------
//     // Scenario 1
//     // ------------------------------------------------------------
//     const MpcScenario scenario1 = {
//         // target 1
//         1.2, 1.8, 0.6, 0.0,

//         // cylinder obstacle
//         -0.6, 2.0, 0.3,

//         // sphere obstacle
//         0.6, 1.8, 0.6, 0.3
//     };

//     // ------------------------------------------------------------
//     // Scenario 2
//     // ------------------------------------------------------------
//     const MpcScenario scenario2 = {
//         // target 2
//         1.0, 2.5, 0.9, 0.0,

//         // cylinder obstacle
//          0.7, 1.3, 0.3,

//         // sphere obstacle
//         0.37, 1.82, 0.3, 0.3
//     };

//     // ------------------------------------------------------------
//     // Scenario 3
//     //
//     // Same target as scenario 2.
//     // No obstacles.
//     //
//     // Obstacles are disabled by:
//     //   radius = 0.0
//     //   center far away
//     //
//     // This avoids adding meaningful obstacle constraints while keeping
//     // the same C API signature.
//     // ------------------------------------------------------------
//     const MpcScenario scenario3 = {
//         // same target as scenario 2
//         scenario2.pxf,
//         scenario2.pyf,
//         scenario2.pzf,
//         scenario2.psif,

//         // disabled cylinder obstacle
//         1000.0, 1000.0, 0.0,

//         // disabled sphere obstacle
//         1000.0, 1000.0, 1000.0, 0.0
//     };

//     // ------------------------------------------------------------
//     // Output CSV
//     // ------------------------------------------------------------
//     const std::string output_csv =
//         "mpc_" + std::to_string(num_mpc_iterations)
//       + "_step_extracted_states.csv";

//     std::ofstream csv(output_csv);

//     if (!csv.is_open()) {
//         std::cerr << "Failed to open " << output_csv << "\n";
//         return 1;
//     }

//     csv << "iteration,"
//         << "scenario,"
//         << "initial_px,initial_py,initial_pz,initial_psi,"
//         << "initial_vx,initial_vy,initial_vz,initial_w,"
//         << "target_x,target_y,target_z,target_psi,"
//         << "cyl_x,cyl_y,cyl_radius,"
//         << "sphere_x,sphere_y,sphere_z,sphere_radius,"
//         << "sample_tau,sample_t,"
//         << "extracted_x,extracted_y,extracted_z,extracted_psi,"
//         << "extracted_vx,extracted_vy,extracted_vz,extracted_w,"
//         << "objective\n";

//     std::cout << std::fixed << std::setprecision(10);
//     csv       << std::fixed << std::setprecision(10);

//     std::cout << "Running receding-horizon OCP for "
//               << num_mpc_iterations << " iterations\n";

//     std::cout << "Scenario 1 iterations: "
//               << scenario1_iterations << "\n";

//     std::cout << "Scenario 2 iterations: "
//               << scenario2_iterations << "\n";

//     std::cout << "Scenario 3 iterations: "
//               << scenario3_iterations << "\n";

//     std::cout << "Scenario 1 -> 2 switch after iteration "
//               << scenario1_iterations << "\n";

//     std::cout << "Scenario 2 -> 3 switch after iteration "
//               << scenario1_iterations + scenario2_iterations << "\n";

//     std::cout << "N = " << N << ", tf = " << tf << "\n";

//     std::cout << "Ts = " << Ts
//               << ", t_sample = " << t_sample
//               << ", tau_sample = " << tau_sample << "\n\n";

//     for (int iter = 0; iter < num_mpc_iterations; ++iter) {
//         const int iteration_number = iter + 1;

//         const int scenario2_start = scenario1_iterations;
//         const int scenario3_start = scenario1_iterations + scenario2_iterations;

//         const bool use_scenario1 = iter < scenario2_start;
//         const bool use_scenario2 = iter >= scenario2_start && iter < scenario3_start;
//         const bool use_scenario3 = iter >= scenario3_start;

//         const MpcScenario& scenario =
//             use_scenario1 ? scenario1 :
//             use_scenario2 ? scenario2 :
//                             scenario3;

//         const int scenario_id =
//             use_scenario1 ? 1 :
//             use_scenario2 ? 2 :
//                             3;

//         if ((iter == scenario2_start || iter == scenario3_start) &&
//             reset_velocity_when_scenario_changes) {
//             vx_cur = 0.0;
//             vy_cur = 0.0;
//             vz_cur = 0.0;
//             w_cur  = 0.0;
//         }

//         const double initial_px  = px_cur;
//         const double initial_py  = py_cur;
//         const double initial_pz  = pz_cur;
//         const double initial_psi = psi_cur;

//         const double initial_vx = vx_cur;
//         const double initial_vy = vy_cur;
//         const double initial_vz = vz_cur;
//         const double initial_w  = w_cur;

//         std::cout << "============================================================\n";
//         std::cout << "MPC iteration " << iteration_number << " / "
//                   << num_mpc_iterations << "\n";

//         std::cout << "Scenario: " << scenario_id << "\n";

//         std::cout << "Initial state: "
//                   << "px=" << initial_px
//                   << ", py=" << initial_py
//                   << ", pz=" << initial_pz
//                   << ", psi=" << initial_psi
//                   << ", vx=" << initial_vx
//                   << ", vy=" << initial_vy
//                   << ", vz=" << initial_vz
//                   << ", w=" << initial_w
//                   << "\n";

//         std::cout << "Target: "
//                   << "x=" << scenario.pxf
//                   << ", y=" << scenario.pyf
//                   << ", z=" << scenario.pzf
//                   << ", psi=" << scenario.psif
//                   << "\n";

//         std::cout << "Cylinder obstacle: "
//                   << "x=" << scenario.cyl_x
//                   << ", y=" << scenario.cyl_y
//                   << ", radius=" << scenario.cyl_radius
//                   << "\n";

//         std::cout << "Sphere obstacle: "
//                   << "x=" << scenario.sphere_x
//                   << ", y=" << scenario.sphere_y
//                   << ", z=" << scenario.sphere_z
//                   << ", radius=" << scenario.sphere_radius
//                   << "\n";

//         // --------------------------------------------------------
//         // Create problem
//         // --------------------------------------------------------
//         PointSetProblem* prob = create_point_set_problem(
//             N, tf,

//             px_max, px_min,
//             py_max, py_min,
//             pz_max, pz_min,
//             psi_max, psi_min,

//             v_max, w_max,
//             a_max, aw_max,

//             px_cur, py_cur, pz_cur, psi_cur,
//             vx_cur, vy_cur, vz_cur, w_cur,

//             scenario.pxf,
//             scenario.pyf,
//             scenario.pzf,
//             scenario.psif,

//             scenario.cyl_x,
//             scenario.cyl_y,
//             scenario.cyl_radius,

//             scenario.sphere_x,
//             scenario.sphere_y,
//             scenario.sphere_z,
//             scenario.sphere_radius
//         );

//         if (!prob) {
//             std::cerr << "create_point_set_problem() returned nullptr at iteration "
//                       << iteration_number << "\n";
//             csv.close();
//             return 1;
//         }

//         // --------------------------------------------------------
//         // Solve
//         // --------------------------------------------------------
//         solve_point_set_problem(prob);

//         const double J = get_final_objective_value(prob);

//         std::vector<double> solution(solution_size, 0.0);
//         get_solution(prob, solution.data(), solution_size);

//         // --------------------------------------------------------
//         // Extract coefficient blocks
//         // --------------------------------------------------------
//         std::vector<double> px_coeffs;
//         std::vector<double> py_coeffs;
//         std::vector<double> pz_coeffs;
//         std::vector<double> psi_coeffs;

//         std::vector<double> vx_coeffs;
//         std::vector<double> vy_coeffs;
//         std::vector<double> vz_coeffs;
//         std::vector<double> w_coeffs;

//         const bool ok_px  = extract_block(solution, 0, Np1, px_coeffs);
//         const bool ok_py  = extract_block(solution, 1, Np1, py_coeffs);
//         const bool ok_pz  = extract_block(solution, 2, Np1, pz_coeffs);
//         const bool ok_psi = extract_block(solution, 3, Np1, psi_coeffs);

//         const bool ok_vx = extract_block(solution, 4, Np1, vx_coeffs);
//         const bool ok_vy = extract_block(solution, 5, Np1, vy_coeffs);
//         const bool ok_vz = extract_block(solution, 6, Np1, vz_coeffs);
//         const bool ok_w  = extract_block(solution, 7, Np1, w_coeffs);

//         if (!ok_px || !ok_py || !ok_pz || !ok_psi ||
//             !ok_vx || !ok_vy || !ok_vz || !ok_w) {
//             std::cerr << "Failed to extract one or more solution blocks at iteration "
//                       << iteration_number << "\n";

//             destroy_point_set_problem(prob);
//             csv.close();
//             return 1;
//         }

//         // --------------------------------------------------------
//         // Optional: print coefficient vectors for debugging
//         // --------------------------------------------------------
//         print_vector("px_coeffs", px_coeffs);
//         print_vector("py_coeffs", py_coeffs);
//         print_vector("pz_coeffs", pz_coeffs);
//         print_vector("psi_coeffs", psi_coeffs);

//         // --------------------------------------------------------
//         // Sample optimized trajectories
//         // --------------------------------------------------------
//         const double x_cmd = bernstein_eval(
//             px_coeffs,
//             tau_sample,
//             -1.0,
//             1.0
//         );

//         const double y_cmd = bernstein_eval(
//             py_coeffs,
//             tau_sample,
//             -1.0,
//             1.0
//         );

//         const double z_cmd = bernstein_eval(
//             pz_coeffs,
//             tau_sample,
//             -1.0,
//             1.0
//         );

//         const double psi_cmd = bernstein_eval(
//             psi_coeffs,
//             tau_sample,
//             -1.0,
//             1.0
//         );

//         const double vx_cmd = bernstein_eval(
//             vx_coeffs,
//             tau_sample,
//             -1.0,
//             1.0
//         );

//         const double vy_cmd = bernstein_eval(
//             vy_coeffs,
//             tau_sample,
//             -1.0,
//             1.0
//         );

//         const double vz_cmd = bernstein_eval(
//             vz_coeffs,
//             tau_sample,
//             -1.0,
//             1.0
//         );

//         const double w_cmd = bernstein_eval(
//             w_coeffs,
//             tau_sample,
//             -1.0,
//             1.0
//         );

//         std::cout << "Extracted sampled state: "
//                   << "x=" << x_cmd
//                   << ", y=" << y_cmd
//                   << ", z=" << z_cmd
//                   << ", psi=" << psi_cmd
//                   << ", vx=" << vx_cmd
//                   << ", vy=" << vy_cmd
//                   << ", vz=" << vz_cmd
//                   << ", w=" << w_cmd
//                   << "\n";

//         std::cout << "Final objective: " << J << "\n";

//         // --------------------------------------------------------
//         // Save row
//         // --------------------------------------------------------
//         csv << iteration_number << ","
//             << scenario_id << ","
//             << initial_px << ","
//             << initial_py << ","
//             << initial_pz << ","
//             << initial_psi << ","
//             << initial_vx << ","
//             << initial_vy << ","
//             << initial_vz << ","
//             << initial_w << ","
//             << scenario.pxf << ","
//             << scenario.pyf << ","
//             << scenario.pzf << ","
//             << scenario.psif << ","
//             << scenario.cyl_x << ","
//             << scenario.cyl_y << ","
//             << scenario.cyl_radius << ","
//             << scenario.sphere_x << ","
//             << scenario.sphere_y << ","
//             << scenario.sphere_z << ","
//             << scenario.sphere_radius << ","
//             << tau_sample << ","
//             << t_sample << ","
//             << x_cmd << ","
//             << y_cmd << ","
//             << z_cmd << ","
//             << psi_cmd << ","
//             << vx_cmd << ","
//             << vy_cmd << ","
//             << vz_cmd << ","
//             << w_cmd << ","
//             << J << "\n";

//         // --------------------------------------------------------
//         // Use extracted state as next OCP initial condition
//         // --------------------------------------------------------
//         px_cur  = x_cmd;
//         py_cur  = y_cmd;
//         pz_cur  = z_cmd;
//         psi_cur = psi_cmd;

//         if (update_velocity_state_from_solution) {
//             vx_cur = vx_cmd;
//             vy_cur = vy_cmd;
//             vz_cur = vz_cmd;
//             w_cur  = w_cmd;
//         } else {
//             vx_cur = 0.0;
//             vy_cur = 0.0;
//             vz_cur = 0.0;
//             w_cur  = 0.0;
//         }

//         destroy_point_set_problem(prob);
//         prob = nullptr;
//     }

//     csv.close();

//     std::cout << "\nDone.\n";
//     std::cout << "Saved extracted MPC states to:\n";
//     std::cout << "  " << output_csv << "\n";

//     return 0;
// }

#include <iostream>
#include <vector>
#include <iomanip>
#include <limits>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <sstream>
#include <unordered_map>
#include <stdexcept>

// ============================================================
// C-API from shared library
// ============================================================
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
        double pxf,    double pyf,    double pzf,    double psif,

        // Cylinder obstacle parameters
        double cyl_x,
        double cyl_y,
        double cyl_radius,

        // Sphere obstacle parameters
        double sphere_x,
        double sphere_y,
        double sphere_z,
        double sphere_radius
    );

    void   solve_point_set_problem(PointSetProblem* problem);
    void   get_solution(PointSetProblem* problem, double* solution, int n);
    double get_final_objective_value(PointSetProblem* problem);
    void   destroy_point_set_problem(PointSetProblem* problem);
}

// ============================================================
// Scenario
// ============================================================
struct MpcScenario {
    double pxf;
    double pyf;
    double pzf;
    double psif;

    double cyl_x;
    double cyl_y;
    double cyl_radius;

    double sphere_x;
    double sphere_y;
    double sphere_z;
    double sphere_radius;
};

// ============================================================
// Input CSV row
// ============================================================
struct TrajectoryPointRow {
    int point_id;
    double arc_length;
    double x;
    double y;
    double z;
};

// ============================================================
// Utility
// ============================================================
static std::string trim(const std::string& s) {
    const auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const auto last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> tokens;
    std::stringstream ss(line);
    std::string item;

    while (std::getline(ss, item, ',')) {
        tokens.push_back(trim(item));
    }

    // Handle trailing comma
    if (!line.empty() && line.back() == ',') {
        tokens.push_back("");
    }

    return tokens;
}

static int find_column_index(
    const std::vector<std::string>& header,
    const std::string& name
) {
    for (int i = 0; i < static_cast<int>(header.size()); ++i) {
        if (trim(header[i]) == name) {
            return i;
        }
    }
    return -1;
}

static std::vector<TrajectoryPointRow> read_input_csv(
    const std::string& input_csv
) {
    std::ifstream in(input_csv);
    if (!in.is_open()) {
        throw std::runtime_error("Failed to open input CSV: " + input_csv);
    }

    std::string line;
    if (!std::getline(in, line)) {
        throw std::runtime_error("Input CSV is empty: " + input_csv);
    }

    const auto header = split_csv_line(line);

    const int idx_point_id   = find_column_index(header, "point_id");
    const int idx_arc_length = find_column_index(header, "arc_length");
    const int idx_x          = find_column_index(header, "x");
    const int idx_y          = find_column_index(header, "y");
    const int idx_z          = find_column_index(header, "z");

    if (idx_point_id < 0 || idx_arc_length < 0 ||
        idx_x < 0 || idx_y < 0 || idx_z < 0) {
        throw std::runtime_error(
            "Input CSV must contain columns: point_id, arc_length, x, y, z"
        );
    }

    std::vector<TrajectoryPointRow> rows;

    while (std::getline(in, line)) {
        if (trim(line).empty()) {
            continue;
        }

        const auto cols = split_csv_line(line);

        const int max_idx = std::max({idx_point_id, idx_arc_length, idx_x, idx_y, idx_z});
        if (static_cast<int>(cols.size()) <= max_idx) {
            throw std::runtime_error("Malformed CSV row: not enough columns.");
        }

        TrajectoryPointRow row{};
        row.point_id   = static_cast<int>(std::llround(std::stod(cols[idx_point_id])));
        row.arc_length = std::stod(cols[idx_arc_length]);
        row.x          = std::stod(cols[idx_x]);
        row.y          = std::stod(cols[idx_y]);
        row.z          = std::stod(cols[idx_z]);

        rows.push_back(row);
    }

    return rows;
}

// ============================================================
// Extract one block from solution
//
// Layout:
//   block 0: px
//   block 1: py
//   block 2: pz
//   block 3: psi
//   block 4: vx
//   block 5: vy
//   block 6: vz
//   block 7: w
//   block 8: ax
//   block 9: ay
//   block 10: az
//   block 11: aw
// ============================================================
static bool extract_block(
    const std::vector<double>& sol,
    int block_idx,
    int Np1,
    std::vector<double>& out
) {
    const int start = block_idx * Np1;
    const int end   = (block_idx + 1) * Np1;

    if (start < 0 || end > static_cast<int>(sol.size()) || end <= start) {
        return false;
    }

    out.assign(sol.begin() + start, sol.begin() + end);
    return true;
}

static void print_vector(
    const std::string& name,
    const std::vector<double>& v
) {
    std::cout << name << " = [";
    for (size_t i = 0; i < v.size(); ++i) {
        std::cout << std::setprecision(10) << v[i];
        if (i + 1 < v.size()) {
            std::cout << ", ";
        }
    }
    std::cout << "]\n";
}

// ============================================================
// Main
// ============================================================
int main() {
    try {
        // ----------------------------------------------------
        // Input / Output files
        // ----------------------------------------------------
        const std::string input_csv  = "smoothed_trajectory_extracted_points_perfect.csv";
        const std::string output_csv = "mpc_vectors_from_smoothed_trajectory.csv";

        // ----------------------------------------------------
        // Read input rows
        // ----------------------------------------------------
        std::vector<TrajectoryPointRow> input_rows = read_input_csv(input_csv);

        if (input_rows.empty()) {
            std::cerr << "No data rows found in input CSV.\n";
            return 1;
        }

        // ----------------------------------------------------
        // Problem size
        // ----------------------------------------------------
        const int    N  = 4;
        const double tf = 0.9;

        const int Np1 = N + 1;
        const int solution_size = 12 * Np1;

        // ----------------------------------------------------
        // Bounds
        // ----------------------------------------------------
        const double px_min = -std::numeric_limits<double>::infinity();
        const double px_max =  std::numeric_limits<double>::infinity();

        const double py_min = -std::numeric_limits<double>::infinity();
        const double py_max =  std::numeric_limits<double>::infinity();

        const double pz_min = 0.0;
        const double pz_max = 5.0;

        const double psi_min = -3.141592653589793;
        const double psi_max =  3.141592653589793;

        const double v_max  = 1.0;
        const double w_max  = 2.0;

        const double a_max  = 1.5;
        const double aw_max = 8.0;

        // ----------------------------------------------------
        // Scenario 1
        // point_id 1 ... 35
        // ----------------------------------------------------
        const MpcScenario scenario1 = {
            // target 1
            1.2, 1.8, 0.6, 0.0,

            // cylinder obstacle
            -0.6, 2.0, 0.2,

            // sphere obstacle
            0.6, 1.8, 0.6, 0.2
        };

        // ----------------------------------------------------
        // Scenario 2
        // point_id 36 ... 70
        // ----------------------------------------------------
        const MpcScenario scenario2 = {
            // target 2
            1.0, 2.5, 0.9, 0.0,

            // cylinder obstacle
            0.6, 1.45, 0.2,

            // sphere obstacle
            0.37, 1.82, 0.6, 0.2
        };

        // ----------------------------------------------------
        // Scenario 3
        // point_id 71 and above
        // Same target as scenario 2, no obstacles
        // ----------------------------------------------------
        const MpcScenario scenario3 = {
            // same target as scenario 2
            scenario2.pxf,
            scenario2.pyf,
            scenario2.pzf,
            scenario2.psif,

            // disabled cylinder obstacle
            1000.0, 1000.0, 0.0,

            // disabled sphere obstacle
            1000.0, 1000.0, 1000.0, 0.0
        };

        // ----------------------------------------------------
        // Open output CSV
        //
        // Columns 1..9:
        //   1 point_id
        //   2 arc_length
        //   3 input_x
        //   4 input_y
        //   5 input_z
        //   6 init_px
        //   7 init_py
        //   8 init_pz
        //   9 init_psi
        //
        // Columns 10...:
        //   px coefficients
        //   py coefficients
        //   pz coefficients
        //   psi coefficients
        //   objective
        //   scenario info
        // ----------------------------------------------------
        std::ofstream csv(output_csv);
        if (!csv.is_open()) {
            std::cerr << "Failed to open output CSV: " << output_csv << "\n";
            return 1;
        }

        csv << std::fixed << std::setprecision(10);
        std::cout << std::fixed << std::setprecision(10);

        // Header
        csv << "point_id,"
            << "arc_length,"
            << "input_x,"
            << "input_y,"
            << "input_z,"
            << "init_px,"
            << "init_py,"
            << "init_pz,"
            << "init_psi,";

        for (int i = 0; i < Np1; ++i) {
            csv << "px_c" << i << ",";
        }
        for (int i = 0; i < Np1; ++i) {
            csv << "py_c" << i << ",";
        }
        for (int i = 0; i < Np1; ++i) {
            csv << "pz_c" << i << ",";
        }
        for (int i = 0; i < Np1; ++i) {
            csv << "psi_c" << i << ",";
        }

        csv << "objective,"
            << "scenario,"
            << "target_x,"
            << "target_y,"
            << "target_z,"
            << "target_psi,"
            << "cyl_x,"
            << "cyl_y,"
            << "cyl_radius,"
            << "sphere_x,"
            << "sphere_y,"
            << "sphere_z,"
            << "sphere_radius\n";

        std::cout << "Read " << input_rows.size()
                  << " trajectory points from " << input_csv << "\n";

        // ----------------------------------------------------
        // Loop over all rows from input CSV
        // ----------------------------------------------------
        for (size_t row_idx = 0; row_idx < input_rows.size(); ++row_idx) {
            const auto& row = input_rows[row_idx];

            // ------------------------------------------------
            // Choose scenario by point_id
            // ------------------------------------------------
            int scenario_id = 0;
            MpcScenario scenario{};

            if (row.point_id < 36) {
                scenario_id = 1;
                scenario = scenario1;
            } else if (row.point_id < 71) {
                scenario_id = 2;
                scenario = scenario2;
            } else {
                scenario_id = 3;
                scenario = scenario3;
            }

            // ------------------------------------------------
            // Initial condition from CSV row
            // psi always 0
            // velocities always 0
            // ------------------------------------------------
            const double px_cur  = row.x;
            const double py_cur  = row.y;
            const double pz_cur  = row.z;
            const double psi_cur = 0.0;

            const double vx_cur = 0.0;
            const double vy_cur = 0.0;
            const double vz_cur = 0.0;
            const double w_cur  = 0.0;

            std::cout << "============================================================\n";
            std::cout << "Row " << (row_idx + 1) << " / " << input_rows.size()
                      << " | point_id = " << row.point_id
                      << " | scenario = " << scenario_id << "\n";

            std::cout << "Initial: "
                      << "px=" << px_cur
                      << ", py=" << py_cur
                      << ", pz=" << pz_cur
                      << ", psi=" << psi_cur << "\n";

            std::cout << "Target: "
                      << "x=" << scenario.pxf
                      << ", y=" << scenario.pyf
                      << ", z=" << scenario.pzf
                      << ", psi=" << scenario.psif << "\n";

            // ------------------------------------------------
            // Create problem
            // ------------------------------------------------
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

                scenario.pxf,
                scenario.pyf,
                scenario.pzf,
                scenario.psif,

                scenario.cyl_x,
                scenario.cyl_y,
                scenario.cyl_radius,

                scenario.sphere_x,
                scenario.sphere_y,
                scenario.sphere_z,
                scenario.sphere_radius
            );

            if (!prob) {
                std::cerr << "create_point_set_problem() returned nullptr for point_id "
                          << row.point_id << "\n";
                csv.close();
                return 1;
            }

            // ------------------------------------------------
            // Solve
            // ------------------------------------------------
            solve_point_set_problem(prob);

            const double J = get_final_objective_value(prob);

            std::vector<double> solution(solution_size, 0.0);
            get_solution(prob, solution.data(), solution_size);

            // ------------------------------------------------
            // Extract coefficient blocks
            // ------------------------------------------------
            std::vector<double> px_coeffs;
            std::vector<double> py_coeffs;
            std::vector<double> pz_coeffs;
            std::vector<double> psi_coeffs;

            const bool ok_px  = extract_block(solution, 0, Np1, px_coeffs);
            const bool ok_py  = extract_block(solution, 1, Np1, py_coeffs);
            const bool ok_pz  = extract_block(solution, 2, Np1, pz_coeffs);
            const bool ok_psi = extract_block(solution, 3, Np1, psi_coeffs);

            if (!ok_px || !ok_py || !ok_pz || !ok_psi) {
                std::cerr << "Failed to extract one or more coefficient blocks for point_id "
                          << row.point_id << "\n";
                destroy_point_set_problem(prob);
                csv.close();
                return 1;
            }

            // Optional debug print
            print_vector("px_coeffs", px_coeffs);
            print_vector("py_coeffs", py_coeffs);
            print_vector("pz_coeffs", pz_coeffs);
            print_vector("psi_coeffs", psi_coeffs);

            std::cout << "Objective: " << J << "\n";

            // ------------------------------------------------
            // Save one output row
            // ------------------------------------------------
            csv << row.point_id << ","
                << row.arc_length << ","
                << row.x << ","
                << row.y << ","
                << row.z << ","
                << px_cur << ","
                << py_cur << ","
                << pz_cur << ","
                << psi_cur << ",";

            for (int i = 0; i < Np1; ++i) {
                csv << px_coeffs[i] << ",";
            }
            for (int i = 0; i < Np1; ++i) {
                csv << py_coeffs[i] << ",";
            }
            for (int i = 0; i < Np1; ++i) {
                csv << pz_coeffs[i] << ",";
            }
            for (int i = 0; i < Np1; ++i) {
                csv << psi_coeffs[i] << ",";
            }

            csv << J << ","
                << scenario_id << ","
                << scenario.pxf << ","
                << scenario.pyf << ","
                << scenario.pzf << ","
                << scenario.psif << ","
                << scenario.cyl_x << ","
                << scenario.cyl_y << ","
                << scenario.cyl_radius << ","
                << scenario.sphere_x << ","
                << scenario.sphere_y << ","
                << scenario.sphere_z << ","
                << scenario.sphere_radius << "\n";

            destroy_point_set_problem(prob);
            prob = nullptr;
        }

        csv.close();

        std::cout << "\nDone.\n";
        std::cout << "Saved all results to:\n";
        std::cout << "  " << output_csv << "\n";

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}