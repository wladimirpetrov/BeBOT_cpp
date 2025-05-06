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
#include <memory>
#include <filesystem> // C++17 and later
#include <sstream>

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
    const std::vector<double>& get_solution() const {
        return raw_solution_;
    }
    // TNLP overrides
    virtual bool get_nlp_info(Index& n, Index& m,
                            Index& nnz_jac_g, Index& nnz_h_lag,
                            IndexStyleEnum& index_style) override
    {
        Index nSeg     = M_*(N_+1);
        n              = 5*nSeg + 1;
        Index nDyn     = 3*nSeg;
        Index nCont    = (M_-1)*5;
        Index nAccCont = (M_-1)*2;        // NEW: accelerate‐continuity
        Index degEl    = 4*N_;
        Index Lobs     = (2*degEl + 1)*M_;
        Index nObs     = static_cast<Index>(p_obs_.size()/2);
        Index nObsC    = nObs * Lobs;
        m               = nDyn + nCont + nAccCont + nObsC;  // UPDATED total constraints

        std::cout << "[DEBUG get_nlp_info] n="<<n<<" m="<<m<<std::endl;
        nnz_jac_g     = n*m;  // dense for FD
        nnz_h_lag     = 0;
        index_style   = TNLP::C_STYLE;
        return true;
    }

    virtual bool get_bounds_info(Index /*n*/, Number* x_l, Number* x_u,
                               Index /*m*/, Number* g_l, Number* g_u) override
    {
        Index nSeg   = M_*(N_+1);
        Index totalX = 5*nSeg + 1;
        // --- variable bounds ---
        for(Index i=0;i<totalX-1;++i){
        x_l[i]=0.0;             // all control points ≥0
        x_u[i]= tf_;            // ≤ final time
        }
        x_l[0]=x_u[0]=x_init_;
        x_l[nSeg-1]=x_u[nSeg-1]=x_final_;
        Index off2 = nSeg;
        x_l[off2]=x_u[off2]=y_init_;
        x_l[off2+nSeg-1]=x_u[off2+nSeg-1]=y_final_;
        Index off3=2*nSeg;
        x_l[off3]=x_u[off3]=headin_;
        x_l[off3+nSeg-1]=x_u[off3+nSeg-1]=headout_;
        Index off4=3*nSeg;
        for(Index i=off4;i<off4+nSeg;++i){
        x_l[i]=-v_max_; x_u[i]=v_max_;
        }
        Index off5=4*nSeg;
        for(Index i=off5;i<off5+nSeg;++i){
        x_l[i]=-omega_max_; x_u[i]=omega_max_;
        }
        // final time
        x_l[totalX-1]=0.0; x_u[totalX-1]= tf_;

        // --- constraint bounds ---
        Index mDyn = 3*nSeg;
        Index idx=0;
        // 1) dynamics (equality)
        for(Index i=0;i<mDyn;++i){
        g_l[idx]=0.0; g_u[idx]=0.0; ++idx;
        }
        // 2) continuity of x1,x2,psi,V,omega (equality)
        for(int seg=0; seg<M_-1; ++seg){
            for(int c=0; c<5; ++c){
                g_l[idx]=0.0; g_u[idx]=0.0; ++idx;
            }
        }
        // 3) continuity of accelerations (NEW) (equality)
        for(int seg=0; seg<M_-1; ++seg){
            for(int c=0; c<2; ++c){
                g_l[idx]=0.0; g_u[idx]=0.0; ++idx;
            }
        }
        // 4) obstacle‐avoidance (ineq)
        Index degEl = 4*N_;
        Index Lobs  = (2*degEl+1)*M_;
        Index nObs  = p_obs_.size()/2;
        Index nObsC = nObs * Lobs;
        for(Index k=0;k<nObsC;++k){
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
        int dim   = N_+1;
        int nSeg  = M_*dim;
        // 1) slice out variables
        std::vector<double> x1(x,      x + nSeg);
        std::vector<double> x2(x+nSeg, x + 2*nSeg);
        std::vector<double> psi(x+2*nSeg, x + 3*nSeg);
        std::vector<double> V  (x+3*nSeg, x + 4*nSeg);
        std::vector<double> om (x+4*nSeg, x + 5*nSeg);
        // 2) rebuild Bézier and differentiation matrix
        tf_ = x[n-1];
        auto tknots = generateTknots();
        piecewiseBebot_ = PiecewiseBeBOT(N_, tknots);
        piecewiseBebot_.calculate();
        std::vector<double> Dm_flat = piecewiseBebot_.getDifferentiationMatrixFlat();

        // 3) compute Dx1, Dx2, Dpsi via MKL
        std::vector<double> Dx1(nSeg), Dx2(nSeg), Dpsi(nSeg);
        for(int seg=0; seg<M_; ++seg){
        int row0 = seg*dim, col0 = seg*dim;
        double* A = Dm_flat.data() + seg*dim*dim;
        cblas_dgemv(CblasRowMajor, CblasTrans, dim, dim,
                    1.0, A, dim,
                    x1.data()+col0, 1, 0.0, Dx1.data()+row0, 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, dim, dim,
                    1.0, A, dim,
                    x2.data()+col0, 1, 0.0, Dx2.data()+row0, 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, dim, dim,
                    1.0, A, dim,
                    psi.data()+col0,1, 0.0, Dpsi.data()+row0,1);
        }

        // 4) dynamics constraints
        for(int i=0;i<nSeg;++i){
            g[i]           = Dx1[i] - V[i]*std::cos(psi[i]);
            g[nSeg+i]      = Dx2[i] - V[i]*std::sin(psi[i]);
            g[2*nSeg+i]    = Dpsi[i] - om[i];
        }

        // 5) continuity of states (x1,x2,psi,V,omega)
        Index idx = 3*nSeg;
        for(int seg=0; seg<M_-1; ++seg){
            int endN = (seg+1)*dim -1, nxt = endN+1;
            g[idx++] = x1[endN]  - x1[nxt];
            g[idx++] = x2[endN]  - x2[nxt];
            g[idx++] = psi[endN] - psi[nxt];
            g[idx++] = V[endN]   - V[nxt];
            g[idx++] = om[endN]  - om[nxt];
        }

        // 6) continuity of accelerations **NEW**
        std::vector<double> Vdot(nSeg);
        std::vector<double> Omegadot(nSeg);

        // --- DEBUG PRINTS: dump tknots, V, om, and first A‐block ---
        // std::cout << "tknots: ";
        // for (double t : tknots) std::cout << t << "  ";
        // std::cout << "\n";

        // std::cout << "Full V vector: ";
        // for (int i = 0; i < nSeg; ++i) std::cout << V[i] << "  ";
        // std::cout << "\n";

        // std::cout << "Full om vector: ";
        // for (int i = 0; i < nSeg; ++i) std::cout << om[i] << "  ";
        // std::cout << "\n";

        // // ── DEBUG: print full Dm_flat block by block ──
        // std::cout << "Full Dm_flat:\n";
        // for(int seg = 0; seg < M_; ++seg) {
        //     std::cout << "-- Segment " << seg << " --\n";
        //     double* Aseg = Dm_flat.data() + seg * dim * dim;
        //     for(int i = 0; i < dim; ++i) {
        //         for(int j = 0; j < dim; ++j) {
        //             std::cout << Aseg[i*dim + j] << " ";
        //         }
        //         std::cout << "\n";
        //     }
        // }
        // std::cout << std::endl;

        // peek at the first segment’s differentiation matrix block A (dim×dim)
        //int dim = N_+1;
        // double* A0 = Dm_flat.data();  // first seg block
        // std::cout << "A block (segment 0):\n";
        // for (int i = 0; i < dim; ++i) {
        //     for (int j = 0; j < dim; ++j) {
        //         std::cout << std::setw(10) << A0[i*dim + j] << " ";
        //     }
        //     std::cout << "\n";
        // }
        // std::cout << std::flush;


        for(int seg=0; seg<M_; ++seg){
        int row0=seg*dim, col0=seg*dim;
        double* A = Dm_flat.data() + seg*dim*dim;
        cblas_dgemv(CblasRowMajor, CblasTrans, dim, dim,
                    1.0, A, dim,
                    V.data()+col0,1, 0.0, Vdot.data()+row0,1);
        cblas_dgemv(CblasRowMajor, CblasTrans, dim, dim,
                    1.0, A, dim,
                    om.data()+col0,1,0.0,Omegadot.data()+row0,1);
        }

        // 4) print Vdot and Omegadot after differentiation
        // std::cout << "Vdot: ";
        // for(int i=0; i<nSeg; ++i) std::cout << Vdot[i] << " ";
        // std::cout << "\n";
        // std::cout << "Omegadot: ";
        // for(int i=0; i<nSeg; ++i) std::cout << Omegadot[i] << " ";
        // std::cout << "\n";

        for(int seg=0; seg<M_-1; ++seg){
            int endN=(seg+1)*dim -1, nxt=endN+1;
            g[idx++] = Vdot[endN]     - Vdot[nxt];
            g[idx++] = Omegadot[endN] - Omegadot[nxt];
        }

        // 7) obstacle avoidance (inequality)
        int degEl = 4*N_;
        int Lobs  = (2*degEl+1)*M_;
        int nObs  = p_obs_.size()/2;
        double sep2 = sep_*sep_;
        idx = 3*nSeg + (M_-1)*5 + (M_-1)*2;
        auto x1_el = PiecewiseDegElevMatrix(x1, M_, N_, degEl);
        auto x2_el = PiecewiseDegElevMatrix(x2, M_, N_, degEl);
        for(int o=0; o<nObs; ++o){
            double xo = p_obs_[o], yo = p_obs_[nObs+o];
            std::vector<double> dx(x1_el.size()), dy(x1_el.size());
            for(size_t i=0;i<dx.size();++i){
                dx[i]=x1_el[i]-xo; dy[i]=x2_el[i]-yo;
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
        raw_solution_.assign(x, x + n);
        
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
        writeToCSV(final_time_, flattened_result_x1, "data/x1.csv");
        writeToCSV(piecewiseBebot_.getNodes(), solution_x1_, "data/x1_controlpoints.csv");

        // X2
        std::vector<std::vector<double>> solution_x2_d(1, std::vector<double>(solution_x2_.begin(), solution_x2_.end()));
        piecewisebernsteinpoly_result_x2_ = PiecewiseBernsteinPoly(solution_x2_d, tknots, final_time_);
        // Flatten bernsteinpoly_psi_result_
        std::vector<double> flattened_result_x2;
        for (const auto& row : piecewisebernsteinpoly_result_x2_) {
            flattened_result_x2.insert(flattened_result_x2.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result_x2, "data/x2.csv");
        writeToCSV(piecewiseBebot_.getNodes(), solution_x2_, "data/x2_controlpoints.csv");

        // psi
        std::vector<std::vector<double>> solution_psi_d(1, std::vector<double>(solution_psi_.begin(), solution_psi_.end()));
        piecewisebernsteinpoly_result_psi_ = PiecewiseBernsteinPoly(solution_psi_d, tknots, final_time_);
        // Flatten bernsteinpoly_psi_result_
        std::vector<double> flattened_result_psi;
        for (const auto& row : piecewisebernsteinpoly_result_psi_) {
            flattened_result_psi.insert(flattened_result_psi.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result_psi, "data/psi.csv");
        writeToCSV(piecewiseBebot_.getNodes(), solution_psi_, "data/psi_controlpoints.csv");

        // v
        std::vector<std::vector<double>> solution_v_d(1, std::vector<double>(solution_v_.begin(), solution_v_.end()));
        piecewisebernsteinpoly_result_v_ = PiecewiseBernsteinPoly(solution_v_d, tknots, final_time_);
        // Flatten bernsteinpoly_v_result_
        std::vector<double> flattened_result_v;
        for (const auto& row : piecewisebernsteinpoly_result_v_) {
            flattened_result_v.insert(flattened_result_v.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result_v, "data/v.csv");
        writeToCSV(piecewiseBebot_.getNodes(), solution_v_, "data/v_controlpoints.csv");

        // om
        std::vector<std::vector<double>> solution_om_d(1, std::vector<double>(solution_om_.begin(), solution_om_.end()));
        piecewisebernsteinpoly_result_om_ = PiecewiseBernsteinPoly(solution_om_d, tknots, final_time_);
        // Flatten bernsteinpoly_om_result_
        std::vector<double> flattened_result_om;
        for (const auto& row : piecewisebernsteinpoly_result_om_) {
            flattened_result_om.insert(flattened_result_om.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result_om, "data/om.csv");
        writeToCSV(piecewiseBebot_.getNodes(), solution_om_, "data/om_controlpoints.csv");

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
        writeToCSV(cont_times, cont_x1,  "data/x1_continuity.csv");
        writeToCSV(cont_times, cont_x2,  "data/x2_continuity.csv");
        writeToCSV(cont_times, cont_psi, "data/psi_continuity.csv");
        writeToCSV(cont_times, cont_v,   "data/v_continuity.csv");
        writeToCSV(cont_times, cont_om,  "data/om_continuity.csv");


        tf_ = x[n-1];
        //auto tknots = generateTknots();
        //piecewiseBebot_ = PiecewiseBeBOT(N_, tknots);
        //piecewiseBebot_.calculate();
        std::vector<double> Dm_flat = piecewiseBebot_.getDifferentiationMatrixFlat();

        int dim = N_+1;
        solution_vdot_.assign(nSeg, 0.0);
        solution_omegadot_.assign(nSeg, 0.0);
        for(int seg=0; seg<M_; ++seg){
            int row0=seg*dim, col0=row0;
            double* A = Dm_flat.data() + seg*dim*dim;
            cblas_dgemv(CblasRowMajor, CblasTrans, dim, dim,
                        1.0, A, dim,
                        solution_v_.data()+col0, 1, 0.0,
                        solution_vdot_.data()+row0, 1);
            cblas_dgemv(CblasRowMajor, CblasTrans, dim, dim,
                        1.0, A, dim,
                        solution_om_.data()+col0, 1, 0.0,
                        solution_omegadot_.data()+row0, 1);
        }
    
        std::vector<std::vector<double>> solution_vdot_d(1, std::vector<double>(solution_vdot_.begin(), solution_vdot_.end()));
        piecewisebernsteinpoly_result_vdot_ = PiecewiseBernsteinPoly(solution_vdot_d, tknots, final_time_);
        // Flatten bernsteinpoly_vdot_result_
        std::vector<double> flattened_result_vdot;
        for (const auto& row : piecewisebernsteinpoly_result_vdot_) {
            flattened_result_vdot.insert(flattened_result_vdot.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result_vdot, "data/vdot.csv");
        writeToCSV(piecewiseBebot_.getNodes(), solution_vdot_, "data/vdot_controlpoints.csv");

        std::vector<std::vector<double>> solution_omegadot_d(1, std::vector<double>(solution_omegadot_.begin(), solution_omegadot_.end()));
        piecewisebernsteinpoly_result_omegadot_ = PiecewiseBernsteinPoly(solution_omegadot_d, tknots, final_time_);
        // Flatten bernsteinpoly_omegadot_result_
        std::vector<double> flattened_result_omegadot;
        for (const auto& row : piecewisebernsteinpoly_result_omegadot_) {
            flattened_result_omegadot.insert(flattened_result_omegadot.end(), row.begin(), row.end());
        }
        writeToCSV(final_time_, flattened_result_omegadot, "data/omegadot.csv");
        writeToCSV(piecewiseBebot_.getNodes(), solution_omegadot_, "data/omegadot_controlpoints.csv");

        std::vector<double> cont_vdot, cont_omegadot;
        cont_vdot.reserve(M_-1);
        cont_omegadot.reserve(M_-1);
        for (int seg = 1; seg < M_; ++seg) {
            int idx = seg*segSize - 1;
            cont_vdot    .push_back(solution_vdot_[idx]);
            cont_omegadot.push_back(solution_omegadot_[idx]);
        }
        writeToCSV(cont_times, cont_vdot,    "data/vdot_continuity.csv");
        writeToCSV(cont_times, cont_omegadot,"data/omegadot_continuity.csv");

        {
        std::ofstream obsFile("data/obstacles.csv");
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
    std::vector<double> raw_solution_;
    PiecewiseBeBOT piecewiseBebot_;
    std::vector<std::vector<double>> piecewisebernsteinpoly_result_x1_;
    std::vector<std::vector<double>> piecewisebernsteinpoly_result_x2_;
    std::vector<std::vector<double>> piecewisebernsteinpoly_result_psi_;
    std::vector<std::vector<double>> piecewisebernsteinpoly_result_v_;
    std::vector<std::vector<double>> piecewisebernsteinpoly_result_om_;
    std::vector<std::vector<double>> piecewisebernsteinpoly_result_vdot_;
    std::vector<std::vector<double>> piecewisebernsteinpoly_result_omegadot_;
    std::vector<double> solution_x1_, solution_x2_, solution_psi_, solution_v_, solution_om_;
    std::vector<double> solution_vdot_;
    std::vector<double> solution_omegadot_;
    Number final_obj_value_;
    std::vector<double> generateTknots() {
        std::vector<double> t(M_ + 1);
        double dt = tf_ / M_;
        for (int i = 0; i <= M_; ++i) t[i] = i * dt;
        return t;
    }
};

extern "C" {

  /**
   * Create a new Dubins‐car obstacle‐avoidance problem instance.
   * Returns a pointer you must later free with destroy_point_set_problem().
   *
   * @param N             Bézier degree
   * @param M             Number of segments
   * @param tf            Nominal final time
   * @param x_init        Initial x
   * @param x_final       Final x
   * @param y_init        Initial y
   * @param y_final       Final y
   * @param heading       Initial heading
   * @param headout       Final heading
   * @param n_obs         Number of obstacles
   * @param sep           Safety separation distance
   * @param v_max         Maximum speed
   * @param omega_max     Maximum turn rate
   * @param p_obs         Array of obstacle x's followed by y's (length = 2*n_obs)
   * @param p_obs_len     Length of p_obs (must equal 2*n_obs)
   */
  PointSetProblem* create_point_set_problem(
      int N, int M, double tf,
      double x_init, double x_final,
      double y_init, double y_final,
      double heading, double headout,
      double n_obs, double sep,
      double v_max, double omega_max,
      const double* p_obs, int p_obs_len
  ) {
      std::vector<double> obs(p_obs, p_obs + p_obs_len);
      return new PointSetProblem(
          N, M, tf,
          x_init, x_final,
          y_init, y_final,
          heading, headout,
          n_obs, sep,
          v_max, omega_max,
          obs
      );
  }

  /**
   * Solve the given problem.  
   * Returns 0 on success (Solve_Succeeded), non‐zero on failure.
   */
  int solve_point_set_problem(PointSetProblem* problem) {
      SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
      app->Options()->SetStringValue("linear_solver",             "ma57");
      app->Options()->SetStringValue("mu_strategy",               "adaptive");
      app->Options()->SetStringValue("gradient_approximation",    "finite-difference-values");
      app->Options()->SetStringValue("jacobian_approximation",    "finite-difference-values");
      app->Options()->SetStringValue("hessian_approximation",     "limited-memory");
      app->Options()->SetIntegerValue("max_iter",                 5000);
      app->Options()->SetNumericValue("tol",                       1e-6);
      app->RethrowNonIpoptException(true);

      if (app->Initialize() != Solve_Succeeded) {
          return -1;
      }
      ApplicationReturnStatus status = app->OptimizeTNLP(
          SmartPtr<TNLP>(problem)
      );
      return static_cast<int>(status);
  }

  /**
   * Retrieve the raw decision vector:
   * [ x1_control_points..., x2_control_points..., psi..., V..., omega..., final_time ]
   *
   * @param problem  The problem pointer returned by create_point_set_problem()
   * @param out      Preallocated array to hold the result
   * @param len      Length of out; must be >= (5*(N+1)*M + 1)
   */
    void get_solution(PointSetProblem* problem, double* out, int len) {
        const auto& sol = problem->get_solution();   // <— call your get_solution()
        int n = std::min(len, (int)sol.size());
        for (int i = 0; i < n; ++i) {
            out[i] = sol[i];
        }
    }

  /**
   * Return the final objective value (i.e. optimized final time).
   */
  double get_final_objective_value(PointSetProblem* problem) {
      return problem->get_final_obj_value();
  }

  /**
   * Destroy a problem instance and free its memory.
   */
  void destroy_point_set_problem(PointSetProblem* problem) {
      delete problem;
  }

}