#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

extern "C" {
    struct PointSetProblem;

    PointSetProblem* create_point_set_problem(
        int N, double tf,
        double px_max, double px_min,
        double py_max, double py_min,
        double pz_max, double pz_min,
        double psi_max, double psi_min,
        double v_max, double w_max,
        double a_max, double aw_max,
        double px_cur, double py_cur, double pz_cur, double psi_cur,
        double vx_cur, double vy_cur, double vz_cur, double w_cur,
        double pxf, double pyf, double pzf, double psif,
        double sphere_x, double sphere_y, double sphere_z,
        double sphere_radius
    );

    int solve_point_set_problem(PointSetProblem* problem);
    int get_solution_size(PointSetProblem* problem);
    int get_solution(PointSetProblem* problem, double* solution, int n);
    double get_final_objective_value(PointSetProblem* problem);
    void destroy_point_set_problem(PointSetProblem* problem);
}

enum class CaseType {
    BLOCKING,
    NEAR_MISS,
    CLEAR
};

struct ObstacleSample {
    double x;
    double y;
    double z;
    double radius;
    double rho;
    double theta;
    CaseType type;
};

static const char* case_type_name(CaseType type) {
    switch (type) {
        case CaseType::BLOCKING: return "blocking";
        case CaseType::NEAR_MISS: return "near_miss";
        case CaseType::CLEAR: return "clear";
    }
    return "unknown";
}

static double uniform_real(
    std::mt19937_64& rng,
    double lower,
    double upper
) {
    std::uniform_real_distribution<double> dist(lower, upper);
    return dist(rng);
}

static bool file_exists(const std::string& filename) {
    std::ifstream f(filename);
    return f.good();
}

// Never overwrite an existing dataset. If foo.csv exists, use foo_001.csv, etc.
static std::string make_non_overwriting_filename(
    const std::string& requested
) {
    if (!file_exists(requested)) {
        return requested;
    }

    const std::size_t dot = requested.find_last_of('.');
    const std::string stem =
        (dot == std::string::npos) ? requested : requested.substr(0, dot);
    const std::string extension =
        (dot == std::string::npos) ? "" : requested.substr(dot);

    for (int i = 1; i < 100000; ++i) {
        std::ostringstream ss;
        ss << stem << "_" << std::setw(3) << std::setfill('0') << i
           << extension;
        if (!file_exists(ss.str())) {
            return ss.str();
        }
    }

    return stem + "_new" + extension;
}

static bool extract_block(
    const std::vector<double>& solution,
    int block_index,
    int block_size,
    std::vector<double>& block
) {
    const int begin = block_index * block_size;
    const int end = begin + block_size;

    if (begin < 0 || end > static_cast<int>(solution.size())) {
        return false;
    }

    block.assign(solution.begin() + begin, solution.begin() + end);
    return true;
}

// Canonical geometry:
//   start  = [0, 0, 0]
//   target = [1, 0, 0]
// The direct start-target line is the x-axis.
//
// We sample the obstacle geometrically using:
//   x     = longitudinal location along the start-target segment
//   rho   = perpendicular distance from the x-axis
//   theta = azimuth around the x-axis
//
//   y = rho*cos(theta)
//   z = rho*sin(theta)
//
// This avoids bias from sampling y and z independently in a box.
static ObstacleSample generate_obstacle(
    std::mt19937_64& rng,
    CaseType type
) {
    constexpr double PI = 3.14159265358979323846;

    // Sphere radius as a fraction of the unit start-target distance.
    const double radius = uniform_real(rng, 0.06, 0.20);

    // Keep the complete sphere away from start and target.
    const double endpoint_margin = 0.05;
    const double x_min = radius + endpoint_margin;
    const double x_max = 1.0 - radius - endpoint_margin;
    const double x = uniform_real(rng, x_min, x_max);

    double rho = 0.0;

    if (type == CaseType::BLOCKING) {
        // rho < radius -> the sphere intersects the straight-line path.
        // Avoid rho = 0 exactly because that case is rotationally symmetric
        // and can have many equivalent optimal avoidance directions.
        rho = uniform_real(rng, 0.05 * radius, 0.90 * radius);
    }
    else if (type == CaseType::NEAR_MISS) {
        // Does not intersect the line, but stays close to it.
        rho = uniform_real(rng, 1.05 * radius, 1.50 * radius);
    }
    else {
        // Clearly non-blocking cases so the dataset also contains
        // trajectories that should remain relatively direct.
        rho = uniform_real(rng, 1.75 * radius, 3.00 * radius);
    }

    const double theta = uniform_real(rng, 0.0, 2.0 * PI);

    ObstacleSample sample;
    sample.x = x;
    sample.y = rho * std::cos(theta);
    sample.z = rho * std::sin(theta);
    sample.radius = radius;
    sample.rho = rho;
    sample.theta = theta;
    sample.type = type;

    return sample;
}

// Exact requested mixture:
//   50% blocking
//   25% near-miss
//   25% clear
// then shuffled so the file is not ordered by category.
static std::vector<CaseType> make_case_schedule(
    std::size_t num_samples,
    std::mt19937_64& rng
) {
    const std::size_t n_blocking = num_samples / 2;
    const std::size_t n_near = num_samples / 4;
    const std::size_t n_clear = num_samples - n_blocking - n_near;

    std::vector<CaseType> schedule;
    schedule.reserve(num_samples);

    schedule.insert(schedule.end(), n_blocking, CaseType::BLOCKING);
    schedule.insert(schedule.end(), n_near, CaseType::NEAR_MISS);
    schedule.insert(schedule.end(), n_clear, CaseType::CLEAR);

    std::shuffle(schedule.begin(), schedule.end(), rng);
    return schedule;
}

static void write_header(std::ofstream& csv, int N) {
    csv
        << "sample_id,attempt_id,case_type,seed,N,tf,v_max,w_max,a_max,aw_max,"
        << "start_x,start_y,start_z,start_psi,"
        << "start_vx,start_vy,start_vz,start_w,"
        << "target_x,target_y,target_z,target_psi,"
        << "sphere_x,sphere_y,sphere_z,sphere_radius,"
        << "sphere_rho,sphere_theta,sphere_line_gap,objective";

    const int L = N + 1;

    for (int i = 0; i < L; ++i) csv << ",px_cp_" << i;
    for (int i = 0; i < L; ++i) csv << ",py_cp_" << i;
    for (int i = 0; i < L; ++i) csv << ",pz_cp_" << i;
    for (int i = 0; i < L; ++i) csv << ",psi_cp_" << i;

    csv << "\n";
}

int main(int argc, char** argv) {
    // ========================================================
    // Run options
    // ========================================================
    // Default:
    //   10,000 SUCCESSFUL OCP solutions
    //   deterministic seed 123456789
    //   safe non-overwriting CSV filename
    //
    // Examples:
    //   ./single_drone_single_obs_data_collection
    //   ./single_drone_single_obs_data_collection 1000
    //   ./single_drone_single_obs_data_collection 50000 42
    //   ./single_drone_single_obs_data_collection 50000 42 my_data.csv
    std::size_t requested_samples = 10000;
    std::uint64_t seed = 123456789ULL;
    std::string requested_output =
        "single_drone_single_obs_training_dataset.csv";

    if (argc >= 2) {
        requested_samples = static_cast<std::size_t>(std::stoull(argv[1]));
    }
    if (argc >= 3) {
        seed = static_cast<std::uint64_t>(std::stoull(argv[2]));
    }
    if (argc >= 4) {
        requested_output = argv[3];
    }

    if (requested_samples == 0) {
        std::cerr << "ERROR: requested_samples must be > 0.\n";
        return 1;
    }

    // ========================================================
    // OCP constants
    // ========================================================
    const int N = 4;
    const int L = N + 1;
    const int expected_solution_size = 12 * L;
    const double tf = 0.9;

    const double INF = std::numeric_limits<double>::infinity();

    // IMPORTANT:
    // In the canonical frame, z is a rotated coordinate, not necessarily
    // physical world altitude. Therefore pz >= 0 would destroy rotational
    // invariance. All three position axes are unbounded here.
    const double px_min = -INF;
    const double px_max =  INF;
    const double py_min = -INF;
    const double py_max =  INF;
    const double pz_min = -INF;
    const double pz_max =  INF;

    const double psi_min = -3.14159265358979323846;
    const double psi_max =  3.14159265358979323846;

    const double v_max = 1.0;
    const double w_max = 2.0;
    const double a_max = 1.5;
    const double aw_max = 8.0;

    // ========================================================
    // Canonical initial condition: ALWAYS [0, 0, 0]
    // ========================================================
    const double px_cur = 0.0;
    const double py_cur = 0.0;
    const double pz_cur = 0.0;
    const double psi_cur = 0.0;

    const double vx_cur = 0.0;
    const double vy_cur = 0.0;
    const double vz_cur = 0.0;
    const double w_cur  = 0.0;

    // ========================================================
    // Canonical target: ALWAYS [1, 0, 0]
    // ========================================================
    const double pxf = 1.0;
    const double pyf = 0.0;
    const double pzf = 0.0;
    const double psif = 0.0;

    // ========================================================
    // RNG and exact case schedule
    // ========================================================
    std::mt19937_64 rng(seed);
    const std::vector<CaseType> schedule =
        make_case_schedule(requested_samples, rng);

    const std::string output_csv =
        make_non_overwriting_filename(requested_output);

    std::ofstream csv(output_csv);
    if (!csv.is_open()) {
        std::cerr << "ERROR: could not open output CSV: "
                  << output_csv << "\n";
        return 2;
    }

    csv << std::fixed << std::setprecision(12);
    write_header(csv, N);

    const std::size_t n_blocking =
        static_cast<std::size_t>(std::count(
            schedule.begin(), schedule.end(), CaseType::BLOCKING));
    const std::size_t n_near =
        static_cast<std::size_t>(std::count(
            schedule.begin(), schedule.end(), CaseType::NEAR_MISS));
    const std::size_t n_clear =
        static_cast<std::size_t>(std::count(
            schedule.begin(), schedule.end(), CaseType::CLEAR));

    std::cout
        << "============================================================\n"
        << "CANONICAL SINGLE-SPHERE OCP DATA COLLECTION\n"
        << "============================================================\n"
        << "Successful samples requested : " << requested_samples << "\n"
        << "Seed                         : " << seed << "\n"
        << "Output                       : " << output_csv << "\n"
        << "Canonical start              : [0, 0, 0]\n"
        << "Canonical target             : [1, 0, 0]\n"
        << "N                            : " << N << "\n"
        << "Control points per block     : " << L << "\n"
        << "Blocking cases               : " << n_blocking << "\n"
        << "Near-miss cases              : " << n_near << "\n"
        << "Clear cases                  : " << n_clear << "\n"
        << "============================================================\n\n";

    std::size_t total_attempts = 0;
    std::size_t failed_solves = 0;
    const std::size_t max_attempts_per_sample = 50;

    // ========================================================
    // Generate exactly requested_samples SUCCESSFUL solves.
    // If a random case fails, retry another obstacle of the same
    // category so the requested category balance is preserved.
    // ========================================================
    for (std::size_t sample_id = 0;
         sample_id < requested_samples;
         ++sample_id) {

        const CaseType type = schedule[sample_id];
        bool saved = false;

        for (std::size_t local_attempt = 0;
             local_attempt < max_attempts_per_sample;
             ++local_attempt) {

            ++total_attempts;

            const ObstacleSample obstacle =
                generate_obstacle(rng, type);

            PointSetProblem* problem = create_point_set_problem(
                N,
                tf,
                px_max, px_min,
                py_max, py_min,
                pz_max, pz_min,
                psi_max, psi_min,
                v_max, w_max,
                a_max, aw_max,
                px_cur, py_cur, pz_cur, psi_cur,
                vx_cur, vy_cur, vz_cur, w_cur,
                pxf, pyf, pzf, psif,
                obstacle.x,
                obstacle.y,
                obstacle.z,
                obstacle.radius
            );

            if (!problem) {
                std::cerr << "ERROR: create_point_set_problem() returned nullptr.\n";
                return 3;
            }

            const int solve_ok = solve_point_set_problem(problem);

            if (!solve_ok) {
                ++failed_solves;
                destroy_point_set_problem(problem);
                problem = nullptr;
                continue;
            }

            const int solution_size = get_solution_size(problem);
            if (solution_size != expected_solution_size) {
                std::cerr
                    << "ERROR: unexpected solution size. Expected "
                    << expected_solution_size
                    << ", received " << solution_size << ".\n";
                destroy_point_set_problem(problem);
                return 4;
            }

            std::vector<double> solution(solution_size, 0.0);
            const int copied = get_solution(
                problem,
                solution.data(),
                solution_size
            );

            if (copied != solution_size) {
                std::cerr << "ERROR: failed to copy complete solution.\n";
                destroy_point_set_problem(problem);
                return 5;
            }

            const double objective =
                get_final_objective_value(problem);

            std::vector<double> px;
            std::vector<double> py;
            std::vector<double> pz;
            std::vector<double> psi;

            const bool blocks_ok =
                extract_block(solution, 0, L, px) &&
                extract_block(solution, 1, L, py) &&
                extract_block(solution, 2, L, pz) &&
                extract_block(solution, 3, L, psi);

            if (!blocks_ok) {
                std::cerr
                    << "ERROR: failed to extract px/py/pz/psi control points.\n";
                destroy_point_set_problem(problem);
                return 6;
            }

            // ====================================================
            // ONE SUCCESSFUL OCP = ONE CSV ROW
            // ====================================================
            csv
                << sample_id << ","
                << total_attempts << ","
                << case_type_name(type) << ","
                << seed << ","
                << N << ","
                << tf << ","
                << v_max << ","
                << w_max << ","
                << a_max << ","
                << aw_max << ","

                << px_cur << ","
                << py_cur << ","
                << pz_cur << ","
                << psi_cur << ","

                << vx_cur << ","
                << vy_cur << ","
                << vz_cur << ","
                << w_cur << ","

                << pxf << ","
                << pyf << ","
                << pzf << ","
                << psif << ","

                << obstacle.x << ","
                << obstacle.y << ","
                << obstacle.z << ","
                << obstacle.radius << ","
                << obstacle.rho << ","
                << obstacle.theta << ","
                << (obstacle.rho - obstacle.radius) << ","
                << objective;

            for (double v : px)  csv << "," << v;
            for (double v : py)  csv << "," << v;
            for (double v : pz)  csv << "," << v;
            for (double v : psi) csv << "," << v;

            csv << "\n";

            // Save completed work during long runs.
            if ((sample_id + 1) % 100 == 0) {
                csv.flush();
            }

            destroy_point_set_problem(problem);
            problem = nullptr;
            saved = true;
            break;
        }

        if (!saved) {
            std::cerr
                << "\nERROR: sample " << sample_id
                << " could not be solved after "
                << max_attempts_per_sample
                << " resamples of the same category.\n"
                << "Stopping rather than silently changing the dataset distribution.\n";
            csv.close();
            return 7;
        }

        const std::size_t completed = sample_id + 1;
        if (completed <= 10 ||
            completed % 100 == 0 ||
            completed == requested_samples) {
            std::cout
                << "Saved " << completed << " / " << requested_samples
                << " successful samples"
                << " | attempts = " << total_attempts
                << " | failed solves = " << failed_solves
                << "\n";
        }
    }

    csv.flush();
    csv.close();

    std::cout
        << "\n============================================================\n"
        << "DATA COLLECTION COMPLETE\n"
        << "============================================================\n"
        << "Successful samples : " << requested_samples << "\n"
        << "Total attempts      : " << total_attempts << "\n"
        << "Failed solves       : " << failed_solves << "\n"
        << "Saved dataset       : " << output_csv << "\n"
        << "============================================================\n";

    return 0;
}