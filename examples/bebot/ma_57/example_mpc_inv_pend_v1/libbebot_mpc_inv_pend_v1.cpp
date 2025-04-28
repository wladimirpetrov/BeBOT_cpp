#include "../../../../Ipopt_ma57_solver/src/Interfaces/IpIpoptApplication.hpp"
#include "../../../../Ipopt_ma57_solver/src/Interfaces/IpTNLP.hpp"
#include "../../../../include/bebot.h"
#include "../../../../include/bernsteinpoly.h"
#include "mkl.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>
#include <iomanip>
#include <memory>
#include <filesystem> // C++17 and later
#include <sstream>

using namespace Ipopt;

// Helper function to check if a file exists
bool fileExists(const std::string& filename) {
    return std::filesystem::exists(filename);
}

// Function to generate a unique filename with a folder
std::string generateUniqueFilename(const std::string& folder, const std::string& baseName, const std::string& extension) {
    int counter = 0;
    std::string uniqueName = folder + "/" + baseName + extension;
    while (fileExists(uniqueName)) {
        ++counter;
        std::ostringstream oss;
        oss << folder << "/" << baseName << "_" << counter << extension;
        uniqueName = oss.str();
    }
    return uniqueName;
}


// Class Definition
class PointSetProblem : public TNLP {
public:
    PointSetProblem(int N, double tf, double theta0, double thetadot0, double thetaf, double thetadotf) 
        : N_(N), tf_(tf), theta0_(theta0), thetadot0_(thetadot0), thetaf_(thetaf), thetadotf_(thetadotf), 
          theta_max_(5.0), theta_min_(-5.0), thetadot_max_(5.0), thetadot_min_(-5.0),
          f_max_(1000.0), f_min_(-1000.0), bebot_(N, tf) {
        initializeMatrices();
        bebot_.calculate();
        std::filesystem::create_directory("data");
    }

    void initializeMatrices() {
        m_ = 0.5;
        L_ = 1.0;
        g_ = 9.81;
        b_ = 0.05;
        I_ = m_ * L_ * L_;

        A_[0][0] = 0;               A_[0][1] = 1;
        A_[1][0] = m_ * g_ * L_ / I_; A_[1][1] = -b_ / I_;

        B_[0][0] = 0;
        B_[1][0] = 1 / I_;
    }

    // Modify the `writeToCSV` function to use the unique filename generator
    void writeToCSV(const std::vector<double>& times, const std::vector<double>& values, const std::string& baseFilename) {
        std::string filename = generateUniqueFilename("data", baseFilename, ".csv");
        std::ofstream outFile(filename);

        if (!outFile.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        outFile << "Time,Value\n"; // CSV headers
        for (size_t i = 0; i < times.size(); ++i) {
            outFile << std::fixed << std::setprecision(6) << times[i] << "," << values[i] << "\n";
        }

        outFile.close();
    }

    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style) {
        n = 3 * (N_ + 1);
        m = 2 * (N_ + 1);
        nnz_jac_g = n * m;
        nnz_h_lag = 0;
        index_style = TNLP::C_STYLE;
        return true;
    }

    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u) {
        // x vector (corresponding to theta,thetadot)
        for (int i = 0; i < 3 * (N_+ 1); i++) {
            x_l[i] = -std::numeric_limits<double>::infinity();
            x_u[i] = std::numeric_limits<double>::infinity();
        }

        // lower and upper bounds for theta values excluding the 1st and the last ones
        for (int i = 1; i < N_; ++i) {
            x_l[i] = theta_min_;
            x_u[i] = theta_max_;
        }
        // lower and upper bounds for the 1st and the last theta values
        x_l[0] = theta0_;
        x_u[0] = theta0_;
        x_l[N_] = thetaf_;
        x_u[N_] = thetaf_;

        // lower and upper bounds for thetadot values excluding the 1st one
        for (int i = N_ + 2; i < 2 * (N_ + 1); ++i) {
            x_l[i] = thetadot_min_;
            x_u[i] = thetadot_max_;
        }

        // lower and upper bounds for the 1st thetadot value
        x_l[N_ + 1] = thetadot0_; 
        x_u[N_ + 1] = thetadot0_;
        x_l[2 * (N_ + 1) - 1] = thetadotf_;
        x_u[2 * (N_ + 1) - 1] = thetadotf_;

        // lower and upper bounds for f values
        for (int i = 2 * (N_ + 1); i < 3 * (N_ + 1); ++i) {
            x_l[i] = f_min_;
            x_u[i] = f_max_;
        }
        
        // Print constraint bounds if needed
        for (Index i = 0; i < 3 * (N_ + 1); ++i) {
            std::cout << "x_l[" << i << "] = " << x_l[i] << ", x_u[" << i << "] = " << x_u[i] << std::endl;
        }

        // g vector (corresponding to theta_dot,thetadouble_dot)
        
        // first 2 (N_ + 1) vectors theta_dot,_dot,theta_dot,q_dot: -inf to +inf 
        for (int i = 0; i < 2 * (N_ + 1); ++i) {
            g_l[i] = 0;
            g_u[i] = 0;
        }

        // Print constraint bounds if needed
        for (Index i = 0; i < 2 * (N_ + 1); ++i) {
            std::cout << "g_l[" << i << "] = " << g_l[i] << ", g_u[" << i << "] = " << g_u[i] << std::endl;
        }

        return true;
    }


    virtual bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda) {
        for (int i = 0; i < 3 * (N_ + 1); ++i) {
            x[i] = 1.0;
        }
        return true;
    }

    // Objective function
    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
        double w1 = 1.0;
        double w2 = 1.0;
        obj_value = 0.0;

        std::vector<double> thetaf_vector(N_ + 1, thetaf_);
        std::vector<double> thetadotf_vector(N_ + 1, thetadotf_);

        // z_vector and theta_vector from x
        std::vector<double> theta_vector(x, x + (N_ + 1)); // Extract the first N+1 elements for theta
        std::vector<double> thetadot_vector(x + (N_ + 1), x + 2 * (N_ + 1)); // Extract the 2nd N+1 elements for thetadot

        // resulting difference between z and zf as well as theta and thetaf
        std::vector<double> theta_diff(N_ + 1);
        std::vector<double> thetadot_diff(N_ + 1);
        /*
        // Print thetaf_vector and thetadotf_vector
        std::cout << "thetaf_vector: ";
        for (const auto& val : thetaf_vector) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        std::cout << "thetadotf_vector: ";
        for (const auto& val : thetadotf_vector) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        // Print theta_vector and thetadot_vector
        std::cout << "theta_vector: ";
        for (const auto& val : theta_vector) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        std::cout << "thetadot_vector: ";
        for (const auto& val : thetadot_vector) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        */
        vdSub(N_ + 1, theta_vector.data(), thetaf_vector.data(), theta_diff.data()); // theta_diff = theta_vector - thetaf_vector
        vdSub(N_ + 1, thetadot_vector.data(), thetadotf_vector.data(), thetadot_diff.data()); // thetadot_diff = thetadot_vector - thetadotf_vector
        /*
        // Print theta_diff and thetadot_diff
        std::cout << "theta_diff: ";
        for (const auto& val : theta_diff) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        std::cout << "thetadot_diff: ";
        for (const auto& val : thetadot_diff) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        */
        // calculating the squared differences
        vdSqr(N_ + 1, theta_diff.data(), theta_diff.data()); // theta_diff = (theta_vector - thetaf_vector)^2
        vdSqr(N_ + 1, thetadot_diff.data(), thetadot_diff.data()); // thetadot_diff = (thetadot_vector - thetadotf_vector)^2
        /*
        // Print theta_diff^2 and thetadot_diff^2
        std::cout << "theta_diff^2: ";
        for (const auto& val : theta_diff) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        std::cout << "thetadot_diff^2: ";
        for (const auto& val : thetadot_diff) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        */
        // Sum the squared differences using MKL's cblas_dasum
        double theta_diff_sum = cblas_dasum(N_ + 1, theta_diff.data(), 1);
        double thetadot_diff_sum = cblas_dasum(N_ + 1, thetadot_diff.data(), 1);

        //std::cout << "theta_diff_sum: " << theta_diff_sum << std::endl;
        //std::cout << "thetadot_diff_sum: " << thetadot_diff_sum << std::endl;

        // Compute the weighted sums
        double theta_term = w1 * theta_diff_sum;
        double thetadot_term = w2 * thetadot_diff_sum;

        // Sum the weighted terms to get the objective value
        obj_value = theta_term + thetadot_term;
        //std::cout << "theta_term: " << theta_term << std::endl;
        //std::cout << "thetadot_term: " << thetadot_term << std::endl;
        //std::cout << "objective value = " << obj_value << std::endl;


        return true;
    }

    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
        // Get the differentiation matrix
        Bebot Bebot(N_, tf_);
        Bebot.calculate();
        const auto& Dm = Bebot.getDifferentiationMatrix();
        // Printing Dm matrix
        //std::cout << "Dm matrix:\n";
        //for (const auto& row : Dm) {
        //    for (const auto& element : row) {
        //        std::cout << element << " ";
        //    }
        //    std::cout << std::endl;
        //}
        // Extracting vectors from x
        std::vector<double> theta_vector(x, x + (N_ + 1));
        std::vector<double> thetadot_vector(x + (N_ + 1), x + 2 * (N_ + 1));
        std::vector<double> f_vector(x + 2 * (N_ + 1), x + 3 * (N_ + 1));
        /*
        std::cout << "theta_vector: ";
        for (const auto& val : theta_vector) {
            std::cout << val << " ";
        }
        std::cout << "thetadot_vector: ";
        for (const auto& val : thetadot_vector) {
            std::cout << val << " ";
        }
        std::cout << "f_vector: ";
        for (const auto& val : f_vector) {
            std::cout << val << " ";
        }
        */
        
        // Calculate dynamics using differentiation matrix Dm
        std::vector<double> x1dot(N_ + 1);
        std::vector<double> x2dot(N_ + 1);

        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, theta_vector.data(), 1, 0.0, x1dot.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, thetadot_vector.data(), 1, 0.0, x2dot.data(), 1);

        //std::cout << "Matrix A_[1][0]: " << A_[1][0] << ", A_[1][1]: " << A_[1][1] << ", B_[1][0]: " << B_[1][0] << std::endl;

        std::vector<double> A_21_theta(N_ + 1);
        std::vector<double> A_22_thetadot(N_ + 1);
        std::vector<double> B_21_f(N_ + 1);

        for (Index i = 0; i < N_ + 1; ++i) {
            A_21_theta[i] = A_[1][0] * theta_vector[i];
            A_22_thetadot[i] = A_[1][1] * thetadot_vector[i];
            B_21_f[i] = B_[1][0] * f_vector[i];

            //std::cout << "A_21_theta[" << i << "]: " << A_21_theta[i] 
            //          << ", A_22_thetadot[" << i << "]: " << A_22_thetadot[i]
            //          << ", B_21_f[" << i << "]: " << B_21_f[i] << std::endl;
        }

        //for (Index i = 0; i < N_ + 1; ++i) {
        //    std::cout << "x1dot[" << i << "]: " << x1dot[i] 
        //              << ", theta_vector[" << i << "]: " << theta_vector[i] 
        //              << ", x2dot[" << i << "]: " << x2dot[i] << std::endl;
        //}

        for (Index i = 0; i < N_ + 1; ++i) {
            g[i] = x1dot[i] - thetadot_vector[i];
            g[N_ + 1 + i] = x2dot[i] - (A_21_theta[i] + A_22_thetadot[i] + B_21_f[i]);
        }

        //std::cout << "Constraint vector g:" << std::endl;
        //for (Index i = 0; i < 2 * (N_ + 1); ++i) {
        //    std::cout << "g[" << i << "]: " << g[i] << std::endl;
        //}

        // Print g vector
        //std::cout << "g vector: ";
        //for (Index i = 0; i < 2 * (N_ + 1); ++i) {
        //    std::cout << g[i] << " ";
        //}
        //std::cout << std::endl;

            return true;
        }

    // Define the Jacobian of /the constraints
    virtual bool eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index* jCol, Number* values) {
        if (values == NULL) {
            // Return the structure of the Jacobian by setting iRow and jCol
            for (Index i = 0; i < m; i++) {
                for (Index j = 0; j < n; j++) {
                    iRow[i * n + j] = i;
                    jCol[i * n + j] = j;
                }
            }
        }
        return true;
    }

    virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) {
        return true;
    }

    virtual void finalize_solution(
        SolverReturn status, 
        Index n,
        const Number* x,
        const Number* z_L,
        const Number* z_U,
        Index m,
        const Number* g,
        const Number* lambda,
        Number obj_value,
        const IpoptData* ip_data,
        IpoptCalculatedQuantities* ip_cq

    ) { // Resize solution vectors
        solution_x_.resize(3 * (N_ + 1));
        solution_u_.resize(2 * (N_ + 1));
        
        // Extract solution_x_ and solution_u_ from x
        for (Index i = 0; i < 2 * (N_ + 1); ++i) {
            solution_x_[i] = x[i];
        }
        for (Index i = 0; i < (N_ + 1); ++i) {
            solution_u_[i] = x[2 * (N_ + 1) + i];
        }
        
        final_obj_value_ = obj_value;         
        bebot_ = Bebot(N_, tf_);
        bebot_.calculate();
        
        final_time_.resize(1000);
        for (int i = 0; i < 1000; ++i) {
            final_time_[i] = i * tf_ / 999.0;
        }
        
        // Create BernsteinPoly results for each variable
        std::vector<double> theta_vector(solution_x_.begin(), solution_x_.begin() + (N_ + 1));
        std::vector<double> thetadot_vector(solution_x_.begin() + (N_ + 1), solution_x_.begin() + 2 * (N_ + 1));

        std::vector<double> f_vector(solution_u_.begin(), solution_u_.begin() + (N_ + 1));
        
        std::vector<std::vector<double>> theta_2d(1, theta_vector);
        std::vector<std::vector<double>> thetadot_2d(1, thetadot_vector);
        std::vector<std::vector<double>> f_2d(1, f_vector);
 
        
        std::vector<std::vector<double>> bernstein_theta = BernsteinPoly(theta_2d, final_time_, 0, tf_);
        std::vector<std::vector<double>> bernstein_thetadot = BernsteinPoly(thetadot_2d, final_time_, 0, tf_);
        std::vector<std::vector<double>> bernstein_f = BernsteinPoly(f_2d, final_time_, 0, tf_);

   
        // Flatten results for saving
        auto flatten = [](const std::vector<std::vector<double>>& input) {
            std::vector<double> output;
            for (const auto& row : input) {
                output.insert(output.end(), row.begin(), row.end());
            }
            return output;
        };
        writeToCSV(final_time_, flatten(bernstein_theta), "theta.csv");
        writeToCSV(bebot_.getNodes(), theta_vector, "theta_controlpoints.csv");
        writeToCSV(final_time_, flatten(bernstein_thetadot), "thetadot.csv");
        writeToCSV(bebot_.getNodes(), thetadot_vector, "thetadot_controlpoints.csv");
        writeToCSV(final_time_, flatten(bernstein_f), "f.csv");
        writeToCSV(bebot_.getNodes(), f_vector, "f_controlpoints.csv");
       
    }

    const std::vector<Number>& get_solution_x() const { return solution_x_; }
    Number get_final_obj_value() const { return final_obj_value_; }

private:
    int N_;
    double tf_, theta0_, thetadot0_, thetaf_, thetadotf_;
    double theta_max_, theta_min_, thetadot_max_, thetadot_min_, f_max_, f_min_;
    double m_, L_, g_, b_, I_;
    double A_[2][2];
    double B_[2][1];
    Bebot bebot_;
    std::vector<Number> solution_u_;
    std::vector<Number> solution_x2_;
    std::vector<Number> solution_x_;
    Number final_obj_value_;
    std::vector<double> final_time_;
    std::vector<std::vector<double>> bernsteinpoly_resultu_;
    std::vector<std::vector<double>> bernsteinpoly_resultx2_;
    std::vector<std::vector<double>> bernsteinpoly_resultz_;
};

// Exported function
extern "C" {
    void solve_point_set_problem(int N, double tf, double theta0, double thetadot0, double thetaf, double thetadotf, std::vector<double>& solution) {
        SmartPtr<TNLP> problem = new PointSetProblem(N, tf, theta0, thetadot0, thetaf, thetadotf);
        SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
        app->Options()->SetStringValue("linear_solver", "ma57");
        app->Options()->SetStringValue("mu_strategy", "adaptive");
        app->Options()->SetStringValue("gradient_approximation", "finite-difference-values");
        app->Options()->SetStringValue("jacobian_approximation", "finite-difference-values");
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");
        app->Options()->SetIntegerValue("max_iter", 5000);
        app->Options()->SetNumericValue("tol", 1e-3);
        //app->Options()->SetIntegerValue("print_level", 0); 
        app->RethrowNonIpoptException(true);
        
        ApplicationReturnStatus status = app->Initialize();
        if (status != Solve_Succeeded) {
            std::cerr << "IPOPT initialization failed!" << std::endl;
            return;
        }

        status = app->OptimizeTNLP(problem);
        
        if (status == Solve_Succeeded || status == Solved_To_Acceptable_Level) {
            // Retrieve the optimal solution and objective value from the problem
            const auto& solution_x = static_cast<PointSetProblem*>(GetRawPtr(problem))->get_solution_x();
            Number final_obj_value = static_cast<PointSetProblem*>(GetRawPtr(problem))->get_final_obj_value();
        
            std::cout << "Optimal Solution (x): ";
            for (Index i = 0; i < solution_x.size(); i++) {
                std::cout << solution_x[i] << " ";
            }
            std::cout << std::endl;
            std::cout << "Optimal Objective Value: " << final_obj_value << std::endl;
            
            // Retrieve the Bernstein polynomial result from the problem
            //const auto& bernsteinpoly_result = static_cast<PointSetProblem*>(GetRawPtr(problem))->get_bernsteinpoly_result();
        } else {
            std::cout << "IPOPT optimization failed with status " << status << std::endl;
        }
    }

    void get_solution(PointSetProblem* problem, double* solution, int n) {
        const std::vector<double>& sol = problem->get_solution_x();
        std::copy(sol.begin(), sol.end(), solution);
    }

    double get_final_objective_value(PointSetProblem* problem) {
        return problem->get_final_obj_value();
    }

    void destroy_point_set_problem(PointSetProblem* problem) {
        delete problem;
        //problem = nullptr;
    }
}


// vpetrov@lnx-me002:/local/vol00/home/vpetrov/dev/optimization/BeBOT_cpp_v2/examples/bebot/ma_57/example_mpc_inv_pend_v1$ g++ -std=c++17 -shared -fPIC -o libbebot_mpc_inv_pend_v1.so ~/dev/optimization/BeBOT_cpp_v2/examples/bebot/ma_57/example_mpc_inv_pend_v1/libbebot_mpc_inv_pend_v1.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bebot.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinpoly.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteindifferentialmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinmatrix_a2b.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/degelevmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/nchoosek_mod.cpp -I~/dev/optimization/BeBOT_cpp_v2/include -I./Ipopt/src/ -L./Ipopt/src/.libs -lipopt -L/opt/intel/oneapi/mkl/latest/lib/intel64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -ldl -lm -lpthread -lstdc++
