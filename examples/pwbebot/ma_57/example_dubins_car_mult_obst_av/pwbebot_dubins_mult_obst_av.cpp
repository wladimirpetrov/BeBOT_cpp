#include "../../../../Ipopt_ma57_solver/src/Interfaces/IpIpoptApplication.hpp"
#include "../../../../Ipopt_ma57_solver/src/Interfaces/IpTNLP.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "../../../../include/piecewisebebot.h"
#include "../../../../include/piecewisedegelevmatrix.h"
#include "../../../../include/piecewisebernsteinproduct.h"
#include "../../../../include/piecewisebernsteinpoly.h"
#include "mkl.h"

using namespace Ipopt;

class PointSetProblem : public TNLP {
public:
    // Constructor
    PointSetProblem(int N, int M, double tf,
                    double x_init, double x_final,
                    double y_init, double y_final,
                    double headin, double headout,
                    double n_obs, double sep,
                    double v_max, double omega_max,
                    const std::vector<double>& p_obs)
        : N_(N), M_(M), tf_(tf),
          x_init_(x_init), x_final_(x_final),
          y_init_(y_init), y_final_(y_final),
          headin_(headin), headout_(headout),
          n_obs_(n_obs), sep_(sep),
          v_max_(v_max), omega_max_(omega_max),
          p_obs_(p_obs),
          piecewiseBebot_(N_, generateTknots()) {
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

    // Public accessor for objective
    Number get_final_obj_value() const { return final_obj_value_; }

    // TNLP overrides
    virtual bool get_nlp_info(Index& n, Index& m,
                              Index& nnz_jac_g, Index& nnz_h_lag,
                              IndexStyleEnum& index_style) override {
        Index nSeg  = M_ * (N_ + 1);
        n = 5 * nSeg + 1;
        Index nDyn  = 3 * nSeg;
        Index nCont = (M_ - 1) * 5;
        Index degEl = 4 * N_;
        Index Lobs  = (2 * degEl + 1) * M_;
        Index nObs  = static_cast<Index>(p_obs_.size() / 2);
        Index nObsC = nObs * Lobs;
        m = nDyn + nCont + nObsC;
        std::cout 
        << "[DEBUG get_nlp_info] n = " << n 
        << ", m = " << m 
        << std::endl;
        nnz_jac_g = n * m;
        nnz_h_lag = 0;
        index_style = TNLP::C_STYLE;
        return true;
    }

    virtual bool get_bounds_info(Index /*n*/, Number* x_l, Number* x_u,
                                 Index /*m*/, Number* g_l, Number* g_u) override {
        Index nSeg   = M_ * (N_ + 1);
        Index totalX = 5 * nSeg + 1;
        // variables
        for (Index i = 0; i < totalX; ++i) {
            x_l[i] = 0.0;
            x_u[i] =  tf_; 
        }
        // boundary conditions
        x_l[0] = x_u[0] = x_init_;
        x_l[nSeg-1] = x_u[nSeg-1] = x_final_;
        Index off2 = nSeg;
        x_l[off2] = x_u[off2] = y_init_;
        x_l[off2+nSeg-1] = x_u[off2+nSeg-1] = y_final_;
        Index off3 = 2 * nSeg;
        x_l[off3] = x_u[off3] = headin_;
        x_l[off3+nSeg-1] = x_u[off3+nSeg-1] = headout_;
        Index off4 = 3 * nSeg;
        for (Index i = off4; i < off4 + nSeg; ++i) { x_l[i] = -v_max_; x_u[i] = v_max_; }
        Index off5 = 4 * nSeg;
        for (Index i = off5; i < off5 + nSeg; ++i) { x_l[i] = -omega_max_; x_u[i] = omega_max_; }
        // constraints
        Index nDyn = 3 * nSeg;
        for (Index i = 0; i < nDyn; ++i) { g_l[i] = 0.0; g_u[i] = 0.0; }
        Index idx = nDyn;
        for (int seg = 0; seg < M_-1; ++seg)
            for (int c = 0; c < 5; ++c) { g_l[idx] = 0.0; g_u[idx] = 0.0; ++idx; }
        Index degEl = 4 * N_;
        Index Lobs  = (2 * degEl + 1) * M_;
        Index nObs  = p_obs_.size() / 2;
        Index nObsC = nObs * Lobs;
        for (Index k = 0; k < nObsC; ++k) {
            g_l[idx + k] = -std::numeric_limits<double>::infinity();
            g_u[idx + k] = 0.0;
        }
        return true;
    }

    virtual bool get_starting_point(
        Index n, bool /*init_x*/, Number* x,
        bool /*init_z*/,   Number* /*z_L*/, Number* /*z_U*/,
        Index /*init_lambda*/, bool /*init_zL*/, Number* /*zl*/) override
    {
        const int dim  = N_ + 1;
        const int nSeg = M_ * dim;

        // 1) build the same tknots0 as in MATLAB
        //    (not strictly needed for x1/x2, but shown for clarity)
        std::vector<double> tknots0(M_+1);
        for(int i=0; i<=M_; ++i)
            tknots0[i] = double(i) * tf_ / M_;

        // 2) x1_init and x2_init:
        //    for each segment j=0..M_-1, linearly interpolate between pinit and pfin
        for(int j = 0; j < M_; ++j){
            // segment endpoints in x
            double x_start = x_init_ + (x_final_ - x_init_) * double(j)   / M_;
            double x_end   = x_init_ + (x_final_ - x_init_) * double(j+1) / M_;
            // segment endpoints in y
            double y_start = y_init_ + (y_final_ - y_init_) * double(j)   / M_;
            double y_end   = y_init_ + (y_final_ - y_init_) * double(j+1) / M_;

            for(int i = 0; i < dim; ++i){
                double alpha = double(i) / N_;  // goes 0..1
                x[j*dim + i]         = x_start + alpha * (x_end - x_start);
                x[nSeg + j*dim + i]  = y_start + alpha * (y_end - y_start);
            }
        }

        // 3) psi_init = linspace(headin, headout, nSeg)
        for(int i = 0; i < nSeg; ++i){
            x[2*nSeg + i] = headin_ 
                        + (headout_ - headin_) * double(i) / double(nSeg - 1);
        }

        // 4) V_init = vmax/2
        for(int i = 0; i < nSeg; ++i){
            x[3*nSeg + i] = v_max_ * 0.5;
        }

        // 5) omega_init = 0
        for(int i = 0; i < nSeg; ++i){
            x[4*nSeg + i] = 0.0;
        }

        // 6) final time T_init = tf_
        x[n-1] = tf_;

        return true;
    }



    virtual bool eval_f(Index n, const Number* x, bool, Number& obj) override {
        obj = x[n-1];
        return true;
    }

    virtual bool eval_g(Index n, const Number* x, bool, Index m, Number* g) override {
        // 1) Extract final time and rebuild knots
        tf_ = x[n-1];
        auto tknots = generateTknots();
        // 2) Reconstruct PiecewiseBeBOT
        piecewiseBebot_ = PiecewiseBeBOT(N_, tknots);
        piecewiseBebot_.calculate();
        // 3) Get flat differentiation matrix
        std::vector<double> Dm_flat = piecewiseBebot_.getDifferentiationMatrixFlat();
        int dim  = N_ + 1;
        int nSeg = M_ * dim;

        // 4) Slice out x1, x2, psi, V, om
        std::vector<double> x1(x,      x + nSeg);
        std::vector<double> x2(x+nSeg, x + 2*nSeg);
        std::vector<double> psi(x+2*nSeg, x + 3*nSeg);
        std::vector<double> V  (x+3*nSeg, x + 4*nSeg);
        std::vector<double> om (x+4*nSeg, x + 5*nSeg);

        // 5) Use MKL to compute Dx1, Dx2, Dpsi via three gemv calls per segment
        std::vector<double> Dx1(nSeg), Dx2(nSeg), Dpsi(nSeg);
        for(int seg = 0; seg < M_; ++seg){
            int row0 = seg * dim;
            int col0 = seg * dim;
            double* A = Dm_flat.data() + seg*dim*dim;
            // Dx1[row0:row0+dim] = A^T * x1[col0:col0+dim]
            cblas_dgemv(CblasRowMajor, CblasTrans,
                        dim, dim,
                        1.0, A, dim,
                        x1.data()+col0, 1,
                        0.0, Dx1.data()+row0, 1);
            // Dx2:
            cblas_dgemv(CblasRowMajor, CblasTrans,
                        dim, dim,
                        1.0, A, dim,
                        x2.data()+col0, 1,
                        0.0, Dx2.data()+row0, 1);
            // Dpsi:
            cblas_dgemv(CblasRowMajor, CblasTrans,
                        dim, dim,
                        1.0, A, dim,
                        psi.data()+col0, 1,
                        0.0, Dpsi.data()+row0, 1);
        }

        // 6) Dynamic constraints
        for(int i=0; i<nSeg; ++i){
            g[i]           = Dx1[i] - V[i]   * std::cos(psi[i]);
            g[nSeg+i]      = Dx2[i] - V[i]   * std::sin(psi[i]);
            g[2*nSeg+i]    = Dpsi[i] - om[i];
        }

        // 7) Continuity
        Index idx = 3*nSeg;
        for(int seg=0; seg<M_-1; ++seg){
            int endN = (seg+1)*dim -1, nxt = endN+1;
            g[idx++] = x1[endN]  - x1[nxt];
            g[idx++] = x2[endN]  - x2[nxt];
            g[idx++] = psi[endN] - psi[nxt];
            g[idx++] = V[endN]   - V[nxt];
            g[idx++] = om[endN]  - om[nxt];
        }

        // 8) Obstacle avoidance (unchanged)
        Index degEl = 4 * N_;
        Index Lobs  = (2 * degEl + 1) * M_;
        Index nObs  = p_obs_.size()/2;
        idx = 3*nSeg + (M_-1)*5;
        auto x1_el = PiecewiseDegElevMatrix(x1, M_, N_, degEl);
        auto x2_el = PiecewiseDegElevMatrix(x2, M_, N_, degEl);
        double sep2 = sep_*sep_;
        for(int o=0; o<nObs; ++o){
            double xo = p_obs_[o], yo = p_obs_[nObs+o];
            std::vector<double> dx(x1_el.size()), dy(x2_el.size());
            for(size_t i=0;i<dx.size();++i){
                dx[i] = x1_el[i]-xo; dy[i] = x2_el[i]-yo;
            }
            auto bx = PiecewiseBernsteinProduct(dx, dx, M_, degEl);
            auto by = PiecewiseBernsteinProduct(dy, dy, M_, degEl);
            for(size_t i=0;i<bx.size();++i){
                g[idx++] = -(bx[i]+by[i]) + sep2;
            }
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
    // finalize_solution lives inside the class
    virtual void finalize_solution(SolverReturn,
                                   Index n,
                                   const Number* x,
                                   const Number*, const Number*,
                                   Index, const Number*, const Number*,
                                   Number obj,
                                   const IpoptData*,
                                   IpoptCalculatedQuantities*) override {
        Index nSeg = M_*(N_+1);
        solution_x1_.assign(x, x+nSeg);
        solution_x2_.assign(x+nSeg, x+2*nSeg);
        solution_psi_.assign(x+2*nSeg, x+3*nSeg);
        solution_v_.assign(x+3*nSeg, x+4*nSeg);
        solution_om_.assign(x+4*nSeg, x+5*nSeg);
        final_obj_value_ = obj;
        tf_ = x[n-1];
        std::vector<double> tknots = generateTknots();
        piecewiseBebot_ = PiecewiseBeBOT(N_, tknots);
        piecewiseBebot_.calculate();
        
        final_time_.resize(1000);
        for (int i = 0; i < 1000; ++i) {
            final_time_[i] = i * obj / 999.0;
            //std::cout << "final_time_[" << i << "] = " << final_time_[i] << std::endl;
        }

        // X1
        std::vector<std::vector<double>> solution_x1_d(1, std::vector<double>(solution_x1_.begin(), solution_x1_.end()));
        piecewisebernsteinpoly_result_x1_ = PiecewiseBernsteinPoly(solution_x1_d, tknots, final_time_);
        // Flatten bernsteinpoly_x1_result_
        std::vector<double> flattened_result_x1;
        for (const auto& row : piecewisebernsteinpoly_result_x1_) {
            flattened_result_x1.insert(flattened_result_x1.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result_x1, "x1.csv");
        writeToCSV(piecewiseBebot_.getNodes(), solution_x1_, "x1_controlpoints.csv");

        // X2
        std::vector<std::vector<double>> solution_x2_d(1, std::vector<double>(solution_x2_.begin(), solution_x2_.end()));
        piecewisebernsteinpoly_result_x2_ = PiecewiseBernsteinPoly(solution_x2_d, tknots, final_time_);
        // Flatten bernsteinpoly_psi_result_
        std::vector<double> flattened_result_x2;
        for (const auto& row : piecewisebernsteinpoly_result_x2_) {
            flattened_result_x2.insert(flattened_result_x2.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result_x2, "x2.csv");
        writeToCSV(piecewiseBebot_.getNodes(), solution_x2_, "x2_controlpoints.csv");

        // psi
        std::vector<std::vector<double>> solution_psi_d(1, std::vector<double>(solution_psi_.begin(), solution_psi_.end()));
        piecewisebernsteinpoly_result_psi_ = PiecewiseBernsteinPoly(solution_psi_d, tknots, final_time_);
        // Flatten bernsteinpoly_psi_result_
        std::vector<double> flattened_result_psi;
        for (const auto& row : piecewisebernsteinpoly_result_psi_) {
            flattened_result_psi.insert(flattened_result_psi.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result_psi, "psi.csv");
        writeToCSV(piecewiseBebot_.getNodes(), solution_psi_, "psi_controlpoints.csv");

        // v
        std::vector<std::vector<double>> solution_v_d(1, std::vector<double>(solution_v_.begin(), solution_v_.end()));
        piecewisebernsteinpoly_result_v_ = PiecewiseBernsteinPoly(solution_v_d, tknots, final_time_);
        // Flatten bernsteinpoly_v_result_
        std::vector<double> flattened_result_v;
        for (const auto& row : piecewisebernsteinpoly_result_v_) {
            flattened_result_v.insert(flattened_result_v.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result_v, "v.csv");
        writeToCSV(piecewiseBebot_.getNodes(), solution_v_, "v_controlpoints.csv");

        // om
        std::vector<std::vector<double>> solution_om_d(1, std::vector<double>(solution_om_.begin(), solution_om_.end()));
        piecewisebernsteinpoly_result_om_ = PiecewiseBernsteinPoly(solution_om_d, tknots, final_time_);
        // Flatten bernsteinpoly_om_result_
        std::vector<double> flattened_result_om;
        for (const auto& row : piecewisebernsteinpoly_result_om_) {
            flattened_result_om.insert(flattened_result_om.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result_om, "om.csv");
        writeToCSV(piecewiseBebot_.getNodes(), solution_om_, "om_controlpoints.csv");

        std::vector<double> cont_times, cont_x1, cont_x2, cont_psi, cont_v, cont_om;
        //auto tknots = generateTknots();  // already have this above
        int segSize = N_ + 1;
        for(int seg=1; seg < M_; ++seg){
            // time at the end of segment `seg`
            cont_times.push_back(tknots[seg]);
            // index of the last control point in segment seg-1
            int idx = seg * segSize - 1;
            cont_x1.push_back(solution_x1_[idx]);
            cont_x2.push_back(solution_x2_[idx]);
            cont_psi.push_back(solution_psi_[idx]);
            cont_v.push_back(solution_v_[idx]);
            cont_om.push_back(solution_om_[idx]);
        }

        // now write them out
        writeToCSV(cont_times, cont_x1,  "x1_continuity.csv");
        writeToCSV(cont_times, cont_x2,  "x2_continuity.csv");
        writeToCSV(cont_times, cont_psi, "psi_continuity.csv");
        writeToCSV(cont_times, cont_v,   "v_continuity.csv");
        writeToCSV(cont_times, cont_om,  "om_continuity.csv");

        {
        std::ofstream obsFile("obstacles.csv");
        obsFile << "x,y,radius\n";
        int nObs = static_cast<int>(p_obs_.size())/2;
        for(int o = 0; o < nObs; ++o){
            double xo = p_obs_[o];
            double yo = p_obs_[nObs + o];
            obsFile << std::fixed << std::setprecision(6)
                    << xo << "," << yo << "," << sep_ << "\n";
        }
        obsFile.close();
    }


    }

private:
    int N_, M_;
    double tf_, x_init_, x_final_, y_init_, y_final_;
    double headin_, headout_, n_obs_, sep_;
    double v_max_, omega_max_;
    std::vector<double> p_obs_;
    std::vector<double> final_time_;
    PiecewiseBeBOT piecewiseBebot_;
    std::vector<std::vector<double>> piecewisebernsteinpoly_result_x1_;
    std::vector<std::vector<double>> piecewisebernsteinpoly_result_x2_;
    std::vector<std::vector<double>> piecewisebernsteinpoly_result_psi_;
    std::vector<std::vector<double>> piecewisebernsteinpoly_result_v_;
    std::vector<std::vector<double>> piecewisebernsteinpoly_result_om_;
    std::vector<double> solution_x1_, solution_x2_, solution_psi_, solution_v_, solution_om_;
    Number final_obj_value_;
    std::vector<double> generateTknots() {
        std::vector<double> t(M_ + 1);
        double dt = tf_ / M_;
        for (int i = 0; i <= M_; ++i) t[i] = i * dt;
        return t;
    }
};

int main() {
    int N = 5;
    int M = 5;
    double tf = 10.0;
    double x_init = 0.0;
    double x_final = 10.0;
    double y_init = 0.0;
    double y_final = 10.0;
    double heading = 0.5236;
    double headout = 0.5236;
    double n_obs = 2.0;
    double sep =  0.5;
    double v_max =  5.0;
    double omega_max = 1.0;
    std::vector<double> p_obs{1.5,4,6,8, 2.0,3.5,6.0,8.5};

    SmartPtr<TNLP> prob = new PointSetProblem(
        N, M,
        tf, x_init, x_final,
        y_init, y_final,
        heading, headout,
        n_obs, sep,
        v_max, omega_max,
        p_obs
    );
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
    app->Options()->SetStringValue("linear_solver","ma57");
    app->Options()->SetStringValue("mu_strategy","adaptive");
    app->Options()->SetStringValue("gradient_approximation","finite-difference-values");
    app->Options()->SetStringValue("jacobian_approximation","finite-difference-values");
    app->Options()->SetStringValue("hessian_approximation","limited-memory");
    app->Options()->SetIntegerValue("max_iter",5000);
    app->Options()->SetNumericValue("tol",1e-6);
    app->RethrowNonIpoptException(true);

    if (app->Initialize() != Solve_Succeeded) {
        std::cerr << "IPOPT initialization failed!" << std::endl;
        return -1;
    }

    ApplicationReturnStatus status = app->OptimizeTNLP(prob);

    if (status == Solve_Succeeded || status == Solved_To_Acceptable_Level) {
        auto solver = static_cast<PointSetProblem*>(GetRawPtr(prob));
        std::cout << "Optimization succeeded. Objective="
                  << solver->get_final_obj_value()
                  << std::endl;
    } else {
        std::cout << "Optimization failed." << std::endl;
    }

    return 0;
}

// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp_v2/examples/pwbebot/ma_57/example_2$ g++ -o pwbebot_example2_v2 ~/dev/optimization/BeBOT_cpp_v2/examples/pwbebot/ma_57/example_2/pwbebot_example2.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/piecewisebebot.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/piecewisebernsteinpoly.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinpoly.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteindifferentialmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinmatrix_a2b.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/degelevmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/nchoosek_mod.cpp -I~/dev/optimization/BeBOT_cpp_v2/include -I./Ipopt/src/ -I/opt/intel/oneapi/mkl/latest/include -L./Ipopt/src/.libs -L/opt/intel/oneapi/mkl/latest/lib/intel64 -lipopt -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -ldl -lm -lpthread -lstdc++
// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp_v2/examples/pwbebot/ma_57/example_2$ export LD_LIBRARY_PATH=/usr/local/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp_v2/examples/pwbebot/ma_57/example_2$ ./pwbebot_example2_v2 
