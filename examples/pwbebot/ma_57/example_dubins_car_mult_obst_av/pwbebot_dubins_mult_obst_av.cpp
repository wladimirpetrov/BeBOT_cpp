#include "../../../../Ipopt_ma57_solver/src/Interfaces/IpIpoptApplication.hpp"
#include "../../../../Ipopt_ma57_solver/src/Interfaces/IpTNLP.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "../../../../include/piecewisebebot.h"
#include "../../../../include/piecewisebernsteinpoly.h"
#include <iomanip>
#include "mkl.h"

using namespace Ipopt;

class PointSetProblem : public Ipopt::TNLP {
public:
    PointSetProblem(int N, int M, double tf, double x_init, double x_final, double y_init, double y_final, double headin, double headout, double n_obs, double sep, double v_max, double omega_max, const std::vector<double>& p_obs)
        : N_(N), M_(M), tf_(tf), x_init_(x_init), x_final_(x_final), y_init_(y_init), y_final_(y_final), headin_(headin), headout_(headout), n_obs_(n_obs), sep_(sep), v_max_(v_max), omega_max_(omega_max), 
          piecewiseBebot_(N, generateTknots()) { // Initialize PiecewiseBeBOT in the initializer list
        piecewiseBebot_.calculate();
    }

    // Write data to CSV
    void writeToCSV(const std::vector<double>& times, const std::vector<double>& values, const std::string& filename) {
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
        // Define the number of variables, constraints, and Jacobian/Hessian non-zero elements.
        n = 5 * M_* (N_ + 1) + 1; 
        
        // total number of collocation constraints
        Index nSeg  = M_*(N_+1);
        Index nDyn  = 3*nSeg;            // dyn1, dyn2, dyn3
        Index nCont = (M_-1)*5;          // x,y,psi,V,omega continuity at each of M–1 gaps
        m = nDyn + nCont;      // new total # of constraints
                
        //m = 3 * M_* (N_ + 1);
        nnz_jac_g = n*m;  
        nnz_h_lag = 0; 
        index_style = TNLP::C_STYLE;
        return true;
    }

    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u) {
        
        Index nSeg   = M_ * (N_ + 1);
        Index totalX = 5 * nSeg + 1;

        // -- x bounds (all free by default)
        for(Index i = 0; i < totalX; ++i) {
            x_l[i] = -std::numeric_limits<double>::infinity();
            x_u[i] =  std::numeric_limits<double>::infinity();
        }
        // x1 init & final
        x_l[0]           = x_u[0]           = x_init_;
        x_l[nSeg - 1]    = x_u[nSeg - 1]    = x_final_;
        // x2 init & final
        Index off2 = nSeg;
        x_l[off2]        = x_u[off2]        = y_init_;
        x_l[off2+nSeg-1] = x_u[off2+nSeg-1] = y_final_;
        // psi init & final
        Index off3 = 2 * nSeg;
        x_l[off3]        = x_u[off3]        = headin_;
        x_l[off3+nSeg-1] = x_u[off3+nSeg-1] = headout_;
        // V bounds
        Index off4 = 3 * nSeg;
        for(Index i = off4; i < off4 + nSeg; ++i) {
            x_l[i] = -v_max_;
            x_u[i] =  v_max_;
        }
        // omega bounds
        Index off5 = 4 * nSeg;
        for(Index i = off5; i < off5 + nSeg; ++i) {
            x_l[i] = -omega_max_;
            x_u[i] =  omega_max_;
        }

        // -- constraint bounds
        // dynamics (dyn1, dyn2, dyn3) == 0
        Index nDyn = 3 * nSeg;
        for(Index i = 0; i < nDyn; ++i) {
            g_l[i] = 0.0;
            g_u[i] = 0.0;
        }
        // continuity == 0
        Index idx = nDyn;
        for(int seg = 0; seg < M_ - 1; ++seg) {
            for(int c = 0; c < 5; ++c) {
                g_l[idx] = 0.0;
                g_u[idx] = 0.0;
                ++idx;
            }
        }

        // Print constraint bounds
        for (Index i = 0; i < M_ * 2 * (N_ + 1) + M_-1; ++i) {
            std::cout << "g_l[" << i << "] = " << g_l[i] << ", g_u[" << i << "] = " << g_u[i] << std::endl;
        }

        return true;
    }

    // initialization of the starting point
    virtual bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda) {
        for (Index i = 0; i < n-1; ++i) {
            x[i] = 1.0;         // your guess for all except the last
            }
            x[n-1] = tf_;

        //for (Index i = 0; i < M_ * (N_ + 1) + 1; ++i) {
        //    std::cout << "x[" << i << "] = " << x[i] << std::endl;
        //}

        return true;
    }

    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
        // Objective function
        obj_value = x[5 * M_ * (N_ + 1)];
        //std::cout << "objective value = " << obj_value << std::endl;
        return true;
    }
    ///*
    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
        
        tf_ = x[M_ * (N_ + 1)];
        std::vector<double> tknots = generateTknots();
        PiecewiseBeBOT piecewiseBebot(N_, tknots);
        piecewiseBebot.calculate();
        const auto& Dm_flat = piecewiseBebot.getDifferentiationMatrixFlat();

        // 6b) Slice into state/control vectors
        Index nSeg = M_ * (N_ + 1);
        std::vector<double> x1_vec   (x, x +   nSeg);
        std::vector<double> x2_vec   (x +   nSeg, x + 2* nSeg);
        std::vector<double> psi_vec  (x + 2* nSeg, x + 3* nSeg);
        std::vector<double> V_vec    (x + 3* nSeg, x + 4* nSeg);
        std::vector<double> omega_vec(x + 4* nSeg, x + 5* nSeg);

        // 6c) Compute derivatives D·x1, D·x2, D·psi using MKL
        std::vector<double> Dx1(nSeg,0.0), Dx2(nSeg,0.0), Dpsi(nSeg,0.0);
        for(int seg = 0; seg < M_; ++seg) {
            int blk = seg*(N_+1)*(N_+1);
            int off = seg*(N_+1);
            cblas_dgemv(CblasRowMajor, CblasNoTrans,
                        N_+1, N_+1,
                        1.0, &Dm_flat[blk], N_+1,
                             &x1_vec[off], 1,
                        0.0, &Dx1[off],    1);
            cblas_dgemv(CblasRowMajor, CblasNoTrans,
                        N_+1, N_+1,
                        1.0, &Dm_flat[blk], N_+1,
                             &x2_vec[off], 1,
                        0.0, &Dx2[off],    1);
            cblas_dgemv(CblasRowMajor, CblasNoTrans,
                        N_+1, N_+1,
                        1.0, &Dm_flat[blk], N_+1,
                             &psi_vec[off], 1,
                        0.0, &Dpsi[off],   1);
        }

        // 6d) Dynamics: dyn1, dyn2, dyn3
        for(Index i = 0; i < nSeg; ++i) {
            // dyn1 = Dx1 - V*cos(psi)
            g[i] = Dx1[i] - V_vec[i] * std::cos(psi_vec[i]);
            // dyn2 = Dx2 - V*sin(psi)
            g[nSeg + i] = Dx2[i] - V_vec[i] * std::sin(psi_vec[i]);
            // dyn3 = Dpsi - omega
            g[2*nSeg + i] = Dpsi[i] - omega_vec[i];
        }

        // 6e) Continuity constraints: end-of-seg - start-of-next
        Index idx = 3 * nSeg;
        for(int seg = 0; seg < M_ - 1; ++seg) {
            int endN  = (seg+1)*(N_+1) - 1;
            int nextN = endN + 1;
            g[idx++] = x1_vec[endN]    - x1_vec[nextN];
            g[idx++] = x2_vec[endN]    - x2_vec[nextN];
            g[idx++] = psi_vec[endN]   - psi_vec[nextN];
            g[idx++] = V_vec[endN]     - V_vec[nextN];
            g[idx++] = omega_vec[endN] - omega_vec[nextN];
        }

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

    // Define the gradient of the objective function
    virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) {
        return true;
    }

    // Method to finalize the solution
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

    ) 
    
    // x1 - x
    { solution_x_.resize(n-1);
        for (Index i = 0; i < n-1; ++i) {
            solution_x_[i] = x[i];
            //std::cout << "solution_x_[" << i << "] = " << solution_x_[i] << std::endl;
        }
        final_obj_value_ = obj_value;
        
        // Updating tf_ with the optimized value of tf (which is x[M_ * (N_ + 1))
        tf_ = x[M_ * (N_ + 1)];
        //std::cout << "tf_ = " << tf_ << std::endl;

        // Recalculate the PiecewiseBeBOT points with the updated final time tf_
        std::vector<double> tknots = generateTknots();
        piecewiseBebot_ = PiecewiseBeBOT(N_, tknots);
        piecewiseBebot_.calculate();

    
        // Calculating final time t using obj_value as final tf for BernsteinPoly library
        final_time_.resize(1000);
        for (int i = 0; i < 1000; ++i) {
            final_time_[i] = i * obj_value / 999.0;
            //std::cout << "final_time_[" << i << "] = " << final_time_[i] << std::endl;
        }

        // Calculating final time t using obj_value as final tf for BernsteinPoly library
        std::vector<std::vector<double>> solution_x_2d(1, std::vector<double>(solution_x_.begin(), solution_x_.end()));
        piecewisebernsteinpoly_result_ = PiecewiseBernsteinPoly(solution_x_2d, tknots, final_time_);
        
        // After calculating final_time_ and bernsteinpoly_result_
        // Flatten bernsteinpoly_result_
        std::vector<double> flattened_result;
        for (const auto& row : piecewisebernsteinpoly_result_) {
            flattened_result.insert(flattened_result.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result, "x.csv");
        writeToCSV(piecewiseBebot_.getNodes(), solution_x_, "x_controlpoints.csv");

        // x2 - g1
        solution_x2_.resize(n-1);
        for (Index i = 0; i < n-1; ++i) {
            solution_x2_[i] = g[i];
            //std::cout << "solution_x2_[" << i << "] = " << solution_x2_[i] << std::endl;
        }
        std::vector<std::vector<double>> solution_x2_2d(1, std::vector<double>(solution_x2_.begin(), solution_x2_.end()));      
        piecewisebernsteinpoly_resultx2_ = PiecewiseBernsteinPoly(solution_x2_2d, tknots, final_time_);
        
        std::vector<double> flattened_result1;
        for (const auto& row : piecewisebernsteinpoly_resultx2_) {
            flattened_result1.insert(flattened_result1.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result1, "x1.csv");
        writeToCSV(piecewiseBebot_.getNodes(), solution_x2_, "x1_controlpoints.csv");

        // u - g2
        solution_u_.resize(n-1);
        for (Index i = 0; i < n-1; ++i) {
            solution_u_[i] = g[M_*(N_ + 1) + i];
            //std::cout << "solution_u_[" << i << "] = " << solution_u_[i] << std::endl;
        }
        std::vector<std::vector<double>> solution_u_2d(1, std::vector<double>(solution_u_.begin(), solution_u_.end()));        
        piecewisebernsteinpoly_resultu_ = PiecewiseBernsteinPoly(solution_u_2d, tknots, final_time_);
        
        //std::cout << "solution_u_2d :" << std::endl;
        //for (const auto& row : piecewisebernsteinpoly_resultu_) {
        //    for (const auto& elem : row) {
        //        std::cout << elem << " ";
        //    }
        //    std::cout << std::endl; 
        //}

        std::vector<double> flattened_result2;
        for (const auto& row : piecewisebernsteinpoly_resultu_) {
            flattened_result2.insert(flattened_result2.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result2, "u.csv");
        writeToCSV(piecewiseBebot_.getNodes(), solution_u_, "u_controlpoints.csv");

        // Save continuity points
        std::vector<double> continuity_times;
        std::vector<double> continuity_values_x;
        std::vector<double> continuity_values_x1;
        std::vector<double> continuity_values_u;

        for (int i = 1; i < M_; ++i) {
            // The time for each continuity point is at the end of each segment
            double time = tknots[i]; // Assuming tknots contain the segment boundaries
            int index = i * (N_ + 1) - 1;

            continuity_times.push_back(time);
            continuity_values_x.push_back(solution_x_[index]);
            continuity_values_x1.push_back(solution_x2_[index]);
            continuity_values_u.push_back(solution_u_[index]);
        }

        writeToCSV(continuity_times, continuity_values_x, "x_continuity.csv");
        writeToCSV(continuity_times, continuity_values_x1, "x1_continuity.csv");
        writeToCSV(continuity_times, continuity_values_u, "u_continuity.csv");

    }
    // Getter for the solution
    const std::vector<Number>& get_solution_x() const { return solution_x_; }

    // Getter for the final objective value
    Number get_final_obj_value() const { return final_obj_value_; }

private:
    int N_;
    int M_;
    double tf_;
    double x_init_;
    double x_final_;
    double y_init_;
    double y_final_;
    double headin_;
    double headout_;
    double n_obs_;
    double sep_;
    double v_max_;
    double omega_max_;
    std::vector<double> p_obs_;
    
    PiecewiseBeBOT piecewiseBebot_; 
    std::vector<Number> solution_u_;
    std::vector<Number> solution_x2_;
    std::vector<Number> solution_x_;
    Number final_obj_value_;
    std::vector<double> final_time_;
    std::vector<std::vector<double>> piecewisebernsteinpoly_resultu_;
    std::vector<std::vector<double>> piecewisebernsteinpoly_resultx2_;
    std::vector<std::vector<double>> piecewisebernsteinpoly_result_;

    // Helper function to generate tknots
    std::vector<double> generateTknots() {
        std::vector<double> tknots;
        //std::cout << "tf_ = " << tf_ << std::endl;
        //std::cout << "M = " << M_ << std::endl;
        double interval = tf_ / M_;
        for (int i = 0; i <= M_; ++i) {
            tknots.push_back(i * interval);
        }
        return tknots;
    }
public: 
};
int main() {
    int N = 4;
    int M = 5;
    double tf = 10;
    double x_init = 0.0;
    double x_final = 10.0;
    double y_init = 0.0;
    double y_final = 10.0;
    double heading = 1.0472;
    double headout = 0.5236;
    double n_obs = 2;
    double sep = 0.5;
    double v_max = 5;
    double omega_max = 1;
    std::vector<double> p_obs{ 1.5, 4.0, 6.0, 8.0, 2.0, 3.5, 6.0, 8.5 };

    SmartPtr<TNLP> pointSetProblem = new PointSetProblem(N, M, tf, x_init, x_final, y_init, y_final, heading, headout, n_obs, sep, v_max, omega_max, p_obs);
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetStringValue("linear_solver", "ma57");
    // A smaller number pivots for sparsity, a larger number pivots for stability
    //app->Options()->SetNumericValue("ma57_pivtol", 1e-8);//0.99 // 1e-8 // between 0 and 1
    // Ipopt may increase pivtol as high as ma27_pivtolmax to get a more accurate solution to the linear system
    //app->Options()->SetNumericValue("ma57_pivtolmax", 0.99);//0.99 // 0.0001 // between 0 and 1
    // The initial integer workspace memory = liw_init_factor * memory required by unfactored system. 
    // Ipopt will increase the workspace size by ma27_meminc_factor if required.
    //app->Options()->SetNumericValue("ma57_liw_init_factor", 5.0); // 5.0 has to be 
    // The initial real workspace memory = la_init_factor * memory required by unfactored system. 
    // Ipopt will increase the workspace size by ma27_meminc_factor if required
    //app->Options()->SetNumericValue("ma57_la_init_factor", 5.0); // 5.0
    // If the integer or real workspace is not large enough, Ipopt will increase its size by this factor.
    //app->Options()->SetNumericValue("ma57_meminc_factor", 5.0); // 5.0

    app->Options()->SetStringValue("mu_strategy", "adaptive");
    
    app->Options()->SetStringValue("gradient_approximation", "finite-difference-values");
    app->Options()->SetStringValue("jacobian_approximation", "finite-difference-values");

    // Set the Hessian approximation method to limited-memory
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");

    // Adjust the maximum number of iterations
    app->Options()->SetIntegerValue("max_iter", 5000); // Change to my desired maximum iterations

    // Adjust the convergence tolerance
    app->Options()->SetNumericValue("tol", 1e-6); // Change to my desired tolerance



    app->RethrowNonIpoptException(true);
    ApplicationReturnStatus status = app->Initialize();
    if (status != Solve_Succeeded) {
        std::cout << "IPOPT initialization failed!" << std::endl;
        return -1;
    }

    status = app->OptimizeTNLP(pointSetProblem);

    if (status == Solve_Succeeded || status == Solved_To_Acceptable_Level) {
        // Process optimization results here
        const auto& solution_x = static_cast<PointSetProblem*>(GetRawPtr(pointSetProblem))->get_solution_x();
        Number final_obj_value = static_cast<PointSetProblem*>(GetRawPtr(pointSetProblem))->get_final_obj_value();
    
        std::cout << "Optimal Solution (x): ";
        for (Index i = 0; i < solution_x.size(); i++) {
            std::cout << solution_x[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "Optimal Objective Value: " << final_obj_value << std::endl;
    } else {
        std::cout << "IPOPT optimization failed with status " << status << std::endl;
    }

    return 0;
}

// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp_v2/examples/pwbebot/ma_57/example_2$ g++ -o pwbebot_example2_v2 ~/dev/optimization/BeBOT_cpp_v2/examples/pwbebot/ma_57/example_2/pwbebot_example2.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/piecewisebebot.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/piecewisebernsteinpoly.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinpoly.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteindifferentialmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinmatrix_a2b.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/degelevmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/nchoosek_mod.cpp -I~/dev/optimization/BeBOT_cpp_v2/include -I./Ipopt/src/ -I/opt/intel/oneapi/mkl/latest/include -L./Ipopt/src/.libs -L/opt/intel/oneapi/mkl/latest/lib/intel64 -lipopt -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -ldl -lm -lpthread -lstdc++
// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp_v2/examples/pwbebot/ma_57/example_2$ export LD_LIBRARY_PATH=/usr/local/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp_v2/examples/pwbebot/ma_57/example_2$ ./pwbebot_example2_v2 
