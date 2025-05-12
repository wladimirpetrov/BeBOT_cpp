#include "../../../../Ipopt_ma57_solver/src/Interfaces/IpIpoptApplication.hpp"
#include "../../../../Ipopt_ma57_solver/src/Interfaces/IpTNLP.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "../../../../include/bebot.h"
#include "../../../../include/bernsteinpoly.h"
#include <iomanip>
#include "mkl.h"
#include "state_space_matrices.h"

using namespace Ipopt;

class PointSetProblem : public Ipopt::TNLP {
public:
    PointSetProblem(int N, double tf, double delta_v_max, double delta_v_min, double delta_s_max, double delta_s_min, double delta_m_max, double delta_m_min, double zmax, double zmin, double wmax, double wmin, double thetamax, double thetamin, double qmax, double qmin, double z0, double w0, double theta0, double q0, double delta_v0, double delta_s0, double delta_m0, double zf, double thetaf,
                    double t0, double tend, double psi0,
                    double a11, double a12, double a13, double a14, double a21, double a22, double a23, double a24, double a31, double a32, double a33, double a34, double a41, double a42, double a43, double a44, 
                    double b11, double b12, double b13, double b21, double b22, double b23, double b31, double b32, double b33, double b41, double b42, double b43)
        : N_(N), tf_(tf), delta_v_max_(delta_v_max), delta_v_min_(delta_v_min), delta_s_max_(delta_s_max), delta_s_min_(delta_s_min), delta_m_max_(delta_m_max), delta_m_min_(delta_m_min), zmax_(zmax), zmin_(zmin), wmax_(wmax), wmin_(wmin), thetamax_(thetamax), thetamin_(thetamin), qmax_(qmax), qmin_(qmin), z0_(z0), w0_(w0), theta0_(theta0), q0_(q0), delta_v0_(delta_v0), delta_s0_(delta_s0), delta_m0_(delta_m0), zf_(zf), thetaf_(thetaf),
          t0_(t0), tend_(tend), psi0_(psi0), bebot_(N, tf_) {
        std::cout << "Creating PointSetProblem instance" << std::endl;

        // Construct the A matrix
        A_ = {{
            {a11, a12, a13, a14},
            {a21, a22, a23, a24},
            {a31, a32, a33, a34},
            {a41, a42, a43, a44}
        }};

        // Construct the B matrix
        B_ = {{
            {b11, b12, b13},
            {b21, b22, b23},
            {b31, b32, b33},
            {b41, b42, b43}
        }};
        
        bebot_.calculate();
    } 

    void writeToCSV(const std::vector<double>& times, const std::vector<double>& values, const std::string& filename) {
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }
        outFile << "Time,Value\n";
        for (size_t i = 0; i < times.size(); ++i) {
            outFile << std::fixed << std::setprecision(6) << times[i] << "," << values[i] << "\n";
        }
        outFile.close();
    }

    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style) {
        n = 7 * (N_ + 1); 
        m = 4 * (N_ + 1);
        nnz_jac_g = n * m;  
        nnz_h_lag = 0; 
        index_style = TNLP::C_STYLE;
        return true;
    }

    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u) {
        // Precompute bounds to avoid redundant calculations
        std::vector<double> x_lower(n, -std::numeric_limits<double>::infinity());
        std::vector<double> x_upper(n, std::numeric_limits<double>::infinity());

        for (int i = 1; i < N_; ++i) {
            x_lower[i] = zmin_;
            x_upper[i] = zmax_;
        }
        x_lower[0] = x_upper[0] = z0_;
        x_lower[N_] = x_upper[N_] = zf_;

        for (int i = N_ + 2; i < 2 * (N_ + 1); ++i) {
            x_lower[i] = wmin_;
            x_upper[i] = wmax_;
        }
        x_lower[N_ + 1] = x_upper[N_ + 1] = w0_;

        for (int i = 2 * (N_ + 1) + 1; i < 3 * (N_ + 1) - 1; ++i) {
            x_lower[i] = thetamin_;
            x_upper[i] = thetamax_;
        }
        x_lower[2 * (N_ + 1)] = x_upper[2 * (N_ + 1)] = theta0_;
        x_lower[3 * (N_ + 1) - 1] = x_upper[3 * (N_ + 1) - 1] = thetaf_;

        for (int i = 3 * (N_ + 1) + 1; i < 4 * (N_ + 1); ++i) {
            x_lower[i] = qmin_;
            x_upper[i] = qmax_;
        }

        for (int i = 4 * (N_ + 1) + 1; i < 5 * (N_ + 1); ++i) {
            x_lower[i] = delta_v_min_;
            x_upper[i] = delta_v_max_;
        }
        x_lower[4 * (N_ + 1)] = x_upper[4 * (N_ + 1)] = delta_v0_;

        for (int i = 5 * (N_ + 1) + 1; i < 6 * (N_ + 1); ++i) {
            x_lower[i] = delta_s_min_;
            x_upper[i] = delta_s_max_;
        }
        x_lower[5 * (N_ + 1)] = x_upper[5 * (N_ + 1)] = delta_s0_;

        for (int i = 6 * (N_ + 1) + 1; i < 7 * (N_ + 1); ++i) {
            x_lower[i] = delta_m_min_;
            x_upper[i] = delta_m_max_;
        }
        x_lower[6 * (N_ + 1)] = x_upper[6 * (N_ + 1)] = delta_m0_;

        std::copy(x_lower.begin(), x_lower.end(), x_l);
        std::copy(x_upper.begin(), x_u);

        std::fill(g_l, g_l + m, 0);
        std::fill(g_u, g_u + m, 0);

        return true;
    }

    virtual bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda) {
        std::fill(x, x + n, 1);
        return true;
    }

    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
        const double w1 = 20.0;
        const double w2 = 20.0;

        obj_value = 0.0;

        std::vector<double> zf_vector(N_ + 1, zf_);
        std::vector<double> thetaf_vector(N_ + 1, thetaf_);

        std::vector<double> z_vector(x, x + (N_ + 1));
        std::vector<double> theta_vector(x + 2 * (N_ + 1), x + 3 * (N_ + 1));

        std::vector<double> z_diff(N_ + 1);
        std::vector<double> theta_diff(N_ + 1);

        vdSub(N_ + 1, z_vector.data(), zf_vector.data(), z_diff.data());
        vdSub(N_ + 1, theta_vector.data(), thetaf_vector.data(), theta_diff.data());

        vdSqr(N_ + 1, z_diff.data(), z_diff.data());
        vdSqr(N_ + 1, theta_diff.data(), theta_diff.data());

        double z_diff_sum = cblas_dasum(N_ + 1, z_diff.data(), 1);
        double theta_diff_sum = cblas_dasum(N_ + 1, theta_diff.data(), 1);

        obj_value = w1 * z_diff_sum + w2 * theta_diff_sum;
        return true;
    }

    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
        Bebot Bebot(N_, tf_);
        Bebot.calculate();
        const auto& Dm = Bebot.getDifferentiationMatrix();

        std::vector<double> z_vector(x, x + (N_ + 1));
        std::vector<double> w_vector(x + (N_ + 1), x + 2 * (N_ + 1));
        std::vector<double> theta_vector(x + 2 * (N_ + 1), x + 3 * (N_ + 1));
        std::vector<double> q_vector(x + 3 * (N_ + 1), x + 4 * (N_ + 1));
        std::vector<double> delta_v_vector(x + 4 * (N_ + 1), x + 5 * (N_ + 1));
        std::vector<double> delta_s_vector(x + 5 * (N_ + 1), x + 6 * (N_ + 1));
        std::vector<double> delta_m_vector(x + 6 * (N_ + 1), x + 7 * (N_ + 1));

        std::vector<double> dyn1(N_ + 1);
        std::vector<double> dyn2(N_ + 1);
        std::vector<double> dyn3(N_ + 1);
        std::vector<double> dyn4(N_ + 1);

        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, z_vector.data(), 1, 0.0, dyn1.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, w_vector.data(), 1, 0.0, dyn2.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, theta_vector.data(), 1, 0.0, dyn3.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, q_vector.data(), 1, 0.0, dyn4.data(), 1);

        std::vector<double> X1_matrix_flat(4 * (N_ + 1));
        std::vector<double> U_matrix_flat(3 * (N_ + 1));

        for (Index i = 0; i < (N_ + 1); ++i) {
            X1_matrix_flat[i] = z_vector[i];
            X1_matrix_flat[(N_ + 1) + i] = w_vector[i];
            X1_matrix_flat[2 * (N_ + 1) + i] = theta_vector[i];
            X1_matrix_flat[3 * (N_ + 1) + i] = q_vector[i];
        }

        for (Index i = 0; i < (N_ + 1); ++i) {
            U_matrix_flat[i] = delta_v_vector[i];
            U_matrix_flat[(N_ + 1) + i] = delta_s_vector[i];
            U_matrix_flat[2 * (N_ + 1) + i] = delta_m_vector[i];
        }

        std::vector<double> X2_matrix_flat(4 * (N_ + 1), 0.0);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4, (N_ + 1), 4, 1.0, &A_[0][0], 4, X1_matrix_flat.data(), (N_ + 1), 0.0, X2_matrix_flat.data(), (N_ + 1));
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4, (N_ + 1), 3, 1.0, &B_[0][0], 3, U_matrix_flat.data(), (N_ + 1), 1.0, X2_matrix_flat.data(), (N_ + 1));

        std::vector<double> z_x2(N_ + 1);
        std::vector<double> w_x2(N_ + 1);
        std::vector<double> theta_x2(N_ + 1);
        std::vector<double> q_x2(N_ + 1);

        for (Index i = 0; i < (N_ + 1); ++i) {
            z_x2[i] = X2_matrix_flat[i];
            w_x2[i] = X2_matrix_flat[(N_ + 1) + i];
            theta_x2[i] = X2_matrix_flat[2 * (N_ + 1) + i];
            q_x2[i] = X2_matrix_flat[3 * (N_ + 1) + i];
        }

        std::vector<double> g1(N_ + 1);
        std::vector<double> g2(N_ + 1);
        std::vector<double> g3(N_ + 1);
        std::vector<double> g4(N_ + 1);

        vdSub(N_ + 1, dyn1.data(), z_x2.data(), g1.data());
        vdSub(N_ + 1, dyn2.data(), w_x2.data(), g2.data());
        vdSub(N_ + 1, dyn3.data(), theta_x2.data(), g3.data());
        vdSub(N_ + 1, dyn4.data(), q_x2.data(), g4.data());

        for (Index i = 0; i < (N_ + 1); ++i) {
            g[i] = g1[i];
            g[(N_ + 1) + i] = g2[i];
            g[2 * (N_ + 1) + i] = g3[i];
            g[3 * (N_ + 1) + i] = g4[i];
        }
        return true;
    }

    virtual bool eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index* jCol, Number* values) {
        if (values == NULL) {
            for (Index i = 0; i < m; i++) {
                for (Index j = 0; j < n; j++) {
                    iRow[i * n + j] = i;
                    jCol[i * n + j] = j;
                }
            }
        }
        return true;
    }

    virtual bool eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lambda, bool new_lambda, Index nele_hess, Index* iRow, Index* jCol, Number* values) {
        // Provide a dummy implementation as the Hessian is not used
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
    ) { 
        std::cout << "Finalizing solution" << std::endl;

        solution_x_.resize(4 * (N_ + 1));
        solution_u_.resize(3 * (N_ + 1));
        
        for (Index i = 0; i < 4 * (N_ + 1); ++i) {
            solution_x_[i] = x[i];
        }
        for (Index i = 0; i < 3 * (N_ + 1); ++i) {
            solution_u_[i] = x[4 * (N_ + 1) + i];
        }
        
        final_obj_value_ = obj_value;         
        bebot_ = Bebot(N_, tf_);
        bebot_.calculate();
        
        final_time_.resize(1000);
        for (int i = 0; i < 1000; ++i) {
            final_time_[i] = i * tf_ / 999.0;
        }

        std::vector<double> z_vector(solution_x_.begin(), solution_x_.begin() + (N_ + 1));
        std::vector<double> w_vector(solution_x_.begin() + (N_ + 1), solution_x_.begin() + 2 * (N_ + 1));
        std::vector<double> theta_vector(solution_x_.begin() + 2 * (N_ + 1), solution_x_.begin() + 3 * (N_ + 1));
        std::vector<double> q_vector(solution_x_.begin() + 3 * (N_ + 1), solution_x_.end());
        std::vector<double> delta_v_vector(solution_u_.begin(), solution_u_.begin() + (N_ + 1));
        std::vector<double> delta_s_vector(solution_u_.begin() + (N_ + 1), solution_u_.begin() + 2 * (N_ + 1));
        std::vector<double> delta_m_vector(solution_u_.begin() + 2 * (N_ + 1), solution_u_.end());

        std::vector<std::vector<double>> z_2d(1, z_vector);
        std::vector<std::vector<double>> w_2d(1, w_vector);
        std::vector<std::vector<double>> theta_2d(1, theta_vector);
        std::vector<std::vector<double>> q_2d(1, q_vector);
        std::vector<std::vector<double>> delta_v_2d(1, delta_v_vector);
        std::vector<std::vector<double>> delta_s_2d(1, delta_s_vector);
        std::vector<std::vector<double>> delta_m_2d(1, delta_m_vector);

        std::vector<std::vector<double>> bernstein_z = BernsteinPoly(z_2d, final_time_, 0, tf_);
        std::vector<std::vector<double>> bernstein_w = BernsteinPoly(w_2d, final_time_, 0, tf_);
        std::vector<std::vector<double>> bernstein_theta = BernsteinPoly(theta_2d, final_time_, 0, tf_);
        std::vector<std::vector<double>> bernstein_q = BernsteinPoly(q_2d, final_time_, 0, tf_);
        std::vector<std::vector<double>> bernstein_delta_v = BernsteinPoly(delta_v_2d, final_time_, 0, tf_);
        std::vector<std::vector<double>> bernstein_delta_s = BernsteinPoly(delta_s_2d, final_time_, 0, tf_);
        std::vector<std::vector<double>> bernstein_delta_m = BernsteinPoly(delta_m_2d, final_time_, 0, tf_);

        auto flatten = [](const std::vector<std::vector<double>>& input) {
            std::vector<double> output;
            for (const auto& row : input) {
                output.insert(output.end(), row.begin(), row.end());
            }
            return output;
        };
        writeToCSV(final_time_, flatten(bernstein_z), "z.csv");
        writeToCSV(bebot_.getNodes(), z_vector, "z_controlpoints.csv");
        writeToCSV(final_time_, flatten(bernstein_w), "w.csv");
        writeToCSV(bebot_.getNodes(), w_vector, "w_controlpoints.csv");
        writeToCSV(final_time_, flatten(bernstein_theta), "theta.csv");
        writeToCSV(bebot_.getNodes(), theta_vector, "theta_controlpoints.csv");
        writeToCSV(final_time_, flatten(bernstein_q), "q.csv");
        writeToCSV(bebot_.getNodes(), q_vector, "q_controlpoints.csv");
        writeToCSV(final_time_, flatten(bernstein_delta_v), "delta_v.csv");
        writeToCSV(bebot_.getNodes(), delta_v_vector, "delta_v_controlpoints.csv");
        writeToCSV(final_time_, flatten(bernstein_delta_s), "delta_s.csv");
        writeToCSV(bebot_.getNodes(), delta_s_vector, "delta_s_controlpoints.csv");
        writeToCSV(final_time_, flatten(bernstein_delta_m), "delta_m.csv");
        writeToCSV(bebot_.getNodes(), delta_m_vector, "delta_m_controlpoints.csv");
    }

    const std::vector<Number>& get_solution_x() const { return solution_x_; }
    Number get_final_obj_value() const { return final_obj_value_; }

private:
    int N_;
    double tf_;
    double delta_v_max_;
    double delta_v_min_;
    double delta_s_max_;
    double delta_s_min_;
    double delta_m_max_;
    double delta_m_min_;
    double zmax_;
    double zmin_;
    double wmax_;  
    double wmin_;
    double thetamax_;
    double thetamin_;
    double qmax_;
    double qmin_;
    double z0_;
    double w0_;
    double theta0_;
    double q0_;
    double delta_v0_;
    double delta_s0_;
    double delta_m0_;
    double zf_;
    double thetaf_;
    double t0_;
    double tend_;
    double psi0_;
    std::array<std::array<double, 4>, 4> A_;
    std::array<std::array<double, 3>, 4> B_;
    Bebot bebot_;
    std::vector<Number> solution_u_;
    std::vector<Number> solution_x2_;
    std::vector<Number> solution_x_;
    Number final_obj_value_;
    std::vector<double> final_time_;
    std::vector<std::vector<double>> bernsteinpoly_resultu_;
    std::vector<std::vector<double>> bernsteinpoly_resultx2_;
    std::vector<std::vector<double>> bernsteinpoly_resultz_;

public:
    const std::vector<std::vector<double>>& get_bernsteinpoly_result() const { 
        return bernsteinpoly_resultz_; }
};

extern "C" void optimizePointSetProblem(int N, double tf, double delta_v_max, double delta_v_min, double delta_s_max, double delta_s_min, double delta_m_max, double delta_m_min, double zmax, double zmin, double wmax, double wmin, double thetamax, double thetamin, double qmax, double qmin, double z0, double w0, double theta0, double q0, double delta_v0, double delta_s0, double delta_m0, double zf, double thetaf,
                                        double t0, double tend, double psi0,
                                        double a11, double a12, double a13, double a14, double a21, double a22, double a23, double a24, double a31, double a32, double a33, double a34, double a41, double a42, double a43, double a44, 
                                        double b11, double b12, double b13, double b21, double b22, double b23, double b31, double b32, double b33, double b41, double b42, double b43) {
    std::cout << "Starting optimization" << std::endl;

    SmartPtr<TNLP> pointSetProblem = new PointSetProblem(N, tf, delta_v_max, delta_v_min, delta_s_max, delta_s_min, delta_m_max, delta_m_min, zmax, zmin, wmax, wmin, thetamax, thetamin, qmax, qmin, z0, w0, theta0, q0, delta_v0, delta_s0, delta_m0, zf, thetaf,
                                                          t0, tend, psi0,
                                                          a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34, a41, a42, a43, a44, 
                                                          b11, b12, b13, b21, b22, b23, b31, b32, b33, b41, b42, b43);
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
    app->Options()->SetStringValue("linear_solver", "ma57");
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("gradient_approximation", "finite-difference-values");
    app->Options()->SetStringValue("jacobian_approximation", "finite-difference-values");
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    app->Options()->SetIntegerValue("max_iter", 5000);
    app->Options()->SetNumericValue("tol", 1e-3);
    app->Options()->SetIntegerValue("print_level", 0); 
    app->RethrowNonIpoptException(true);
    ApplicationReturnStatus status = app->Initialize();
    if (status != Solve_Succeeded) {
        std::cout << "IPOPT initialization failed!" << std::endl;
        return;
    }
    status = app->OptimizeTNLP(pointSetProblem);
    if (status == Solve_Succeeded || status == Solved_To_Acceptable_Level) {
        std::cout << "Optimization succeeded" << std::endl;

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
}

extern "C" void get_solution(PointSetProblem* problem, double* solution, int n) {
    const std::vector<double>& sol = problem->get_solution_x();
    std::copy(sol.begin(), sol.end(), solution);
}

extern "C" double get_final_objective_value(PointSetProblem* problem) {
    return problem->get_final_obj_value();
}

extern "C" void destroy_problem(PointSetProblem* problem) {
    delete problem;
}
