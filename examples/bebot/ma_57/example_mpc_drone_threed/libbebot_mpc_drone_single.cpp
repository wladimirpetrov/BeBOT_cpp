#include "../../../../Ipopt_ma57_solver/src/Interfaces/IpIpoptApplication.hpp"
#include "../../../../Ipopt_ma57_solver/src/Interfaces/IpTNLP.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "../../../../include/bebot.h"
#include "../../../../include/bernsteinpoly.h"
#include <iomanip>
#include "mkl.h"

#include <array>
#include <string>
#include <limits>
#include <algorithm>

#include <sstream>

using namespace Ipopt;

class PointSetProblem : public Ipopt::TNLP {
public:
    PointSetProblem(int N, double tf, 
        double px_max, double px_min, double py_max, double py_min, double pz_max, double pz_min, double psi_max, double psi_min, 
        double v_max, double w_max,
        double a_max, double aw_max, 
        double px_cur, double py_cur, double pz_cur, double psi_cur, 
        double vx_cur, double vy_cur, double vz_cur, double w_cur, 
        double pxf, double pyf, double pzf, double psif)
        : N_(N), tf_(tf), 
        px_max_(px_max), px_min_(px_min), py_max_(py_max), py_min_(py_min), pz_max_(pz_max), pz_min_(pz_min), psi_max_(psi_max), psi_min_(psi_min),
        v_max_(v_max), w_max_(w_max),
        a_max_(a_max), aw_max_(aw_max), 
        px_cur_(px_cur), py_cur_(py_cur), pz_cur_(pz_cur), psi_cur_(psi_cur),
        vx_cur_(vx_cur), vy_cur_(vy_cur), vz_cur_(vz_cur), w_cur_(w_cur),
        pxf_(pxf), pyf_(pyf), pzf_(pzf), psif_(psif), bebot_(N, tf_) {
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
        n = 12 * (N_ + 1); // p(x,y,z), psi, v(x,y,z), w, a(x,y,z), aw, 
        m = 10 * (N_ + 1); // 4*(N+1) resp + 4*(N+1) resv + 2*(N+1) cspeed+aspeed 
        nnz_jac_g = n * m;  
        nnz_h_lag = 0; 
        index_style = TNLP::C_STYLE;
        return true;
    }

    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u) {
        // Precompute bounds to avoid redundant calculations
        std::vector<double> x_lower(n, -std::numeric_limits<double>::infinity());
        std::vector<double> x_upper(n, std::numeric_limits<double>::infinity());
        // px
        for (int i = 1; i < N_ + 1; ++i) {
            x_lower[i] = px_min_;
            x_upper[i] = px_max_;
        }
        x_lower[0] = x_upper[0] = px_cur_;
        //x_lower[N_] = x_upper[N_] = zf_;

        // py
        for (int i = N_ + 1 + 1; i < 2 * (N_ + 1); ++i) {
            x_lower[i] = py_min_;
            x_upper[i] = py_max_;
        }
        x_lower[(N_ + 1)] = x_upper[(N_ + 1)] = py_cur_;

        // pz
        for (int i = 2 * (N_ + 1) + 1; i < 3 * (N_ + 1); ++i) {
            x_lower[i] = pz_min_;
            x_upper[i] = pz_max_;
        }
        x_lower[2 * (N_ + 1)] = x_upper[2 * (N_ + 1)] = pz_cur_;

        // psi
        for (int i = 3 * (N_ + 1) + 1; i < 4 * (N_ + 1); ++i) {
            x_lower[i] = psi_min_;
            x_upper[i] = psi_max_;
        }
        x_lower[3 * (N_ + 1)] = x_upper[3 * (N_ + 1)] = psi_cur_;//

        // vx
        for (int i = 4 * (N_ + 1) + 1; i < 5 * (N_ + 1); ++i) {
            x_lower[i] = -v_max_;
            x_upper[i] = v_max_;
        }
        x_lower[4 * (N_ + 1)] = x_upper[4 * (N_ + 1)] = vx_cur_;

        // vy
        for (int i = 5 * (N_ + 1) + 1; i < 6 * (N_ + 1); ++i) {
            x_lower[i] = -v_max_;
            x_upper[i] = v_max_;
        }
        x_lower[5 * (N_ + 1)] = x_upper[5 * (N_ + 1)] = vy_cur_;

        // vz
        for (int i = 6 * (N_ + 1) + 1; i < 7 * (N_ + 1); ++i) {
            x_lower[i] = -v_max_;
            x_upper[i] = v_max_;
        }
        x_lower[6 * (N_ + 1)] = x_upper[6 * (N_ + 1)] = vz_cur_;

        // w
        for (int i = 7 * (N_ + 1) + 1; i < 8 * (N_ + 1); ++i) {
            x_lower[i] = -w_max_;
            x_upper[i] = w_max_;
        }
        x_lower[7 * (N_ + 1)] = x_upper[7 * (N_ + 1)] = w_cur_;


        // control input
        // ax
        for (int i = 8 * (N_ + 1); i < 9 * (N_ + 1); ++i) { 
            x_lower[i] = -a_max_;
            x_upper[i] = a_max_;
        }
        
        // ay
        for (int i = 9 * (N_ + 1); i < 10 * (N_ + 1); ++i) { 
            x_lower[i] = -a_max_;
            x_upper[i] = a_max_;
        }

        // ay
        for (int i = 10 * (N_ + 1); i < 11 * (N_ + 1); ++i) {
            x_lower[i] = -a_max_;
            x_upper[i] = a_max_;
        }

        // aw
        for (int i = 11 * (N_ + 1); i < 12 * (N_ + 1); ++i) {
            x_lower[i] = -aw_max_;
            x_upper[i] = aw_max_;
        }

        std::copy(x_lower.begin(), x_lower.end(), x_l);
        std::copy(x_upper.begin(), x_upper.end(), x_u);

        // --------------------
        // g bounds (constraints)
        // m must be 10*(N_+1)
        // --------------------

        // 0..8*(N+1)-1 : dynamics equalities
        for (int i = 0; i < 8 * (N_ + 1); ++i) {
            g_l[i] = 0.0;
            g_u[i] = 0.0;
        }

        // 8*(N+1)..9*(N+1)-1 : c_speed <= 0
        for (int i = 8 * (N_ + 1); i < 9 * (N_ + 1); ++i) {
            g_l[i] = -std::numeric_limits<double>::infinity();
            g_u[i] = 0.0;
        }

        // 9*(N+1)..10*(N+1)-1 : c_accel <= 0
        for (int i = 9 * (N_ + 1); i < 10 * (N_ + 1); ++i) {
            g_l[i] = -std::numeric_limits<double>::infinity();
            g_u[i] = 0.0;
        }

        for (Index i = 0; i < n; ++i) std::cout << "x_l["<<i<<"]="<<x_l[i]<<", x_u["<<i<<"]="<<x_u[i]<<"\n";
        for (Index i = 0; i < m; ++i) std::cout << "g_l["<<i<<"]="<<g_l[i]<<", g_u["<<i<<"]="<<g_u[i]<<"\n";


        return true;
    }

    virtual bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda) {
        for (Index i = 0; i < n; ++i) {
            x[i] = 1.0;
        }
        // for (Index i = 0; i < n; ++i) std::cout << "x0["<<i<<"]="<<x[i]<<"\n";

        return true;
    }


    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
        

        // ---- OPTIONAL: throttle prints (eval_f is called a lot) ----
        // static int f_call = 0;
        // const bool do_print = (f_call < 1000);   // print first 3 calls only
        // ++f_call;

        // auto print_vec = [&](const std::string& name, const std::vector<double>& v, int max_elems = 10) {
        //     if (!do_print) return;
        //     std::cout << std::fixed << std::setprecision(6);
        //     std::cout << name << " (size=" << v.size() << ") = [";
        //     const int M = std::min<int>(static_cast<int>(v.size()), max_elems);
        //     for (int i = 0; i < M; ++i) {
        //         std::cout << v[i];
        //         if (i + 1 < M) std::cout << ", ";
        //     }
        //     if (static_cast<int>(v.size()) > max_elems) std::cout << ", ...";
        //     std::cout << "]\n";
        // };

        // auto print_scalar = [&](const std::string& name, double v) {
        //     if (!do_print) return;
        //     std::cout << std::fixed << std::setprecision(12);
        //     std::cout << name << " = " << v << "\n";
        // };

        const double w_p = 1;
        const double w_psi = 1;
        const double w_a = 0.01;
        const double w_aw = 0.001;

        obj_value = 0.0;
        // referenced xyz and psi
        std::vector<double> pxf_vector(N_ + 1, pxf_);
        std::vector<double> pyf_vector(N_ + 1, pyf_);
        std::vector<double> pzf_vector(N_ + 1, pzf_);
        std::vector<double> psif_vector(N_ + 1, psif_);
        
        // px,py,pz,psi vectors from x
        std::vector<double> px_vector(x, x + (N_ + 1));
        std::vector<double> py_vector(x + 1 * (N_ + 1), x + 2 * (N_ + 1));
        std::vector<double> pz_vector(x + 2 * (N_ + 1), x + 3 * (N_ + 1));
        std::vector<double> psi_vector(x + 3 * (N_ + 1), x + 4 * (N_ + 1));

        // ---- PRINT: references and current decision vectors ----
        // print_vec("pxf_vector",  pxf_vector);
        // print_vec("pyf_vector",  pyf_vector);
        // print_vec("pzf_vector",  pzf_vector);
        // print_vec("psif_vector", psif_vector);

        // print_vec("px_vector",   px_vector);
        // print_vec("py_vector",   py_vector);
        // print_vec("pz_vector",   pz_vector);
        // print_vec("psi_vector",  psi_vector);

        // vx,vy,vz,w vectors from x
        // std::vector<double> vx_vector(x + 4 * (N_ + 1), x + 5 * (N_ + 1));
        // std::vector<double> vy_vector(x + 5 * (N_ + 1), x + 6 * (N_ + 1));
        // std::vector<double> vz_vector(x + 6 * (N_ + 1), x + 7 * (N_ + 1));
        // std::vector<double> w_vector(x + 7 * (N_ + 1), x + 8 * (N_ + 1));

        // ax,ay,az,aw vectors from x
        std::vector<double> ax_vector(x + 8 * (N_ + 1), x + 9 * (N_ + 1));
        std::vector<double> ay_vector(x + 9 * (N_ + 1), x + 10 * (N_ + 1));
        std::vector<double> az_vector(x + 10 * (N_ + 1), x + 11 * (N_ + 1));
        std::vector<double> aw_vector(x + 11 * (N_ + 1), x + 12 * (N_ + 1));

        // print_vec("ax_vector",  ax_vector);
        // print_vec("ay_vector",  ay_vector);
        // print_vec("az_vector",  az_vector);
        // print_vec("aw_vector", aw_vector);

        // --- errors ---
        std::vector<double> px_diff(N_ + 1);
        std::vector<double> py_diff(N_ + 1);
        std::vector<double> pz_diff(N_ + 1);
        std::vector<double> psi_diff(N_ + 1);

        vdSub(N_ + 1, px_vector.data(),     pxf_vector.data(),     px_diff.data());
        vdSub(N_ + 1, py_vector.data(),     pyf_vector.data(),     py_diff.data());
        vdSub(N_ + 1, pz_vector.data(),     pzf_vector.data(),     pz_diff.data());
        vdSub(N_ + 1, psi_vector.data(),     psif_vector.data(),     psi_diff.data());


        // ---- PRINT: diffs after vdSub ----
        // print_vec("px_diff",   px_diff);
        // print_vec("py_diff",   py_diff);
        // print_vec("pz_diff",   pz_diff);
        // print_vec("psi_diff",  psi_diff);

        // --- squared errors ---
        std::vector<double> px_diff_sqr(N_ + 1);
        std::vector<double> py_diff_sqr(N_ + 1);
        std::vector<double> pz_diff_sqr(N_ + 1);
        std::vector<double> psi_diff_sqr(N_ + 1);



        vdSqr(N_ + 1, px_diff.data(),     px_diff_sqr.data());
        vdSqr(N_ + 1, py_diff.data(),     py_diff_sqr.data());
        vdSqr(N_ + 1, pz_diff.data(),     pz_diff_sqr.data());
        vdSqr(N_ + 1, psi_diff.data(),     psi_diff_sqr.data());

        // ---- PRINT: squared diffs ----
        // print_vec("px_diff_sqr",  px_diff_sqr);
        // print_vec("py_diff_sqr",  py_diff_sqr);
        // print_vec("pz_diff_sqr",  pz_diff_sqr);
        // print_vec("psi_diff_sqr", psi_diff_sqr);

        // --- r_tau (MATLAB else-branch) ---
        // tau = linspace(-1, 0.999, N+1)
        // r_tau = 2 ./ (1 - tau).^2
        std::vector<double> r_tau(N_ + 1);
        if (N_ + 1 == 1) {
            const double tau0 = -1.0;
            const double denom = (1.0 - tau0);
            r_tau[0] = 2.0 / (denom * denom);
        } else {
            const double tau0 = -1.0;
            const double tau1 = 0.9;
            const double step = (tau1 - tau0) / static_cast<double>(N_); // (N+1)-1 = N
            for (int i = 0; i < (N_ + 1); ++i) {
                const double tau = tau0 + step * static_cast<double>(i);
                const double denom = (1.0 - tau);
                r_tau[i] = 2.0 / (denom * denom);
            }
        }
        // print_vec("r_tau", r_tau);
        // for (int i = 0; i < N_ + 1; ++i) std::cout << "r_tau[" << i << "]=" << r_tau[i] << "\n";

        // --- elementwise weighted squared errors: r_tau .* (error.^2) ---
        std::vector<double> px_weighted(N_ + 1);
        std::vector<double> py_weighted(N_ + 1);
        std::vector<double> pz_weighted(N_ + 1);
        std::vector<double> psi_weighted(N_ + 1);


        vdMul(N_ + 1, r_tau.data(), px_diff_sqr.data(),     px_weighted.data());
        vdMul(N_ + 1, r_tau.data(), py_diff_sqr.data(),     py_weighted.data());
        vdMul(N_ + 1, r_tau.data(), pz_diff_sqr.data(),     pz_weighted.data());
        vdMul(N_ + 1, r_tau.data(), psi_diff_sqr.data(),     psi_weighted.data());

        // print_vec("px_weighted",  px_weighted);
        // print_vec("py_weighted",  py_weighted);
        // print_vec("pz_weighted",  pz_weighted);
        // print_vec("psi_weighted", psi_weighted);

        // --- sum() via dot with ones ---
        std::vector<double> ones(N_ + 1, 1.0);

        const double sum_px_state     = cblas_ddot(N_ + 1, px_weighted.data(),     1, ones.data(), 1);
        const double sum_py_state     = cblas_ddot(N_ + 1, py_weighted.data(),     1, ones.data(), 1);
        const double sum_pz_state     = cblas_ddot(N_ + 1, pz_weighted.data(),     1, ones.data(), 1);
        const double sum_psi_state     = cblas_ddot(N_ + 1, psi_weighted.data(),     1, ones.data(), 1);

        // ---- PRINT: sums ----
        // print_scalar("sum_px_state",  sum_px_state);
        // print_scalar("sum_py_state",  sum_py_state);
        // print_scalar("sum_pz_state",  sum_pz_state);
        // print_scalar("sum_psi_state", sum_psi_state);

        const double state_term =
            w_p     * sum_px_state + w_p * sum_py_state + w_p * sum_pz_state + w_psi * sum_psi_state;
        // ---- PRINT: state term ----
        // print_scalar("state_term", state_term);
        // --- control term: sum(u.^2) (UNWEIGHTED, same as MATLAB snippet) ---

        std::vector<double> ax_sqr(N_ + 1);
        std::vector<double> ay_sqr(N_ + 1);
        std::vector<double> az_sqr(N_ + 1);
        std::vector<double> aw_sqr(N_ + 1);

        vdSqr(N_ + 1, ax_vector.data(), ax_sqr.data());
        vdSqr(N_ + 1, ay_vector.data(), ay_sqr.data());
        vdSqr(N_ + 1, az_vector.data(), az_sqr.data());
        vdSqr(N_ + 1, aw_vector.data(), aw_sqr.data());

        // ---- PRINT: control squares ----
        // print_vec("ax_sqr", ax_sqr);
        // print_vec("ay_sqr", ay_sqr);
        // print_vec("az_sqr", az_sqr);
        // print_vec("aw_sqr", aw_sqr);


        const double sum_ax = cblas_ddot(N_ + 1, ax_sqr.data(), 1, ones.data(), 1);
        const double sum_ay = cblas_ddot(N_ + 1, ay_sqr.data(), 1, ones.data(), 1);
        const double sum_az = cblas_ddot(N_ + 1, az_sqr.data(), 1, ones.data(), 1);
        const double sum_aw = cblas_ddot(N_ + 1, aw_sqr.data(), 1, ones.data(), 1);

        // ---- PRINT: control sums ----
        // print_scalar("sum_ax", sum_ax);
        // print_scalar("sum_ay", sum_ay);
        // print_scalar("sum_az", sum_az);
        // print_scalar("sum_aw", sum_aw);

        const double control_term = w_a * sum_ax + w_a * sum_ay + w_a * sum_az + w_aw * sum_aw;
        // print_scalar("control_term", control_term);
        obj_value = state_term + control_term;

        return true;
    }

    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {

        // ---- OPTIONAL: throttle prints (eval_g is called a lot) ----
        static int g_call = 0;
        const bool do_print = (g_call < 1000);   // print first 2 calls only
        ++g_call;

        // auto print_vec = [&](const std::string& name, const std::vector<double>& v, int max_elems = 10) {
        //     if (!do_print) return;
        //     std::cout << std::fixed << std::setprecision(6);
        //     std::cout << name << " (size=" << v.size() << ") = [";
        //     const int M = std::min<int>(static_cast<int>(v.size()), max_elems);
        //     for (int i = 0; i < M; ++i) {
        //         std::cout << v[i];
        //         if (i + 1 < M) std::cout << ", ";
        //     }
        //     if (static_cast<int>(v.size()) > max_elems) std::cout << ", ...";
        //     std::cout << "]\n";
        // };

        // auto print_scalar = [&](const std::string& name, double v) {
        //     if (!do_print) return;
        //     std::cout << std::fixed << std::setprecision(12);
        //     std::cout << name << " = " << v << "\n";
        // };

        // auto print_mat = [&](const std::string& name,
        //                     const std::vector<double>& A,
        //                     int rows, int cols,
        //                     int max_rows = 8, int max_cols = 8) {
        //     if (!do_print) return;
        //     std::cout << std::fixed << std::setprecision(6);
        //     std::cout << name << " (" << rows << "x" << cols << ")\n";
        //     const int R = std::min(rows, max_rows);
        //     const int C = std::min(cols, max_cols);

        //     // NOTE: you use Dm with CblasColMajor later.
        //     // So print as column-major: A[r + rows*c]
        //     for (int r = 0; r < R; ++r) {
        //         for (int c = 0; c < C; ++c) {
        //             std::cout << A[r + rows * c] << "  ";
        //         }
        //         if (cols > max_cols) std::cout << "...";
        //         std::cout << "\n";
        //     }
        //     if (rows > max_rows) std::cout << "...\n";
        // };

        Bebot Bebot(N_, tf_);
        Bebot.calculate();
        const auto& Dm = Bebot.getDifferentiationMatrix();

        // ---- PRINT: Dm (partial by default) ----
        // print_mat("Dm", Dm, N_ + 1, N_ + 1);

        // -----------------------------
        // 1) Build tau and r_tau_row
        // -----------------------------
        std::vector<double> tau_vec(N_ + 1);
        std::vector<double> r_tau_row(N_ + 1);

        if (N_ + 1 == 1) {
            tau_vec[0] = -1.0;
            const double denom = 1.0 - tau_vec[0];
            r_tau_row[0] = 2.0 / (denom * denom);
        } else {
            const double tau0 = -1.0;
            const double tau1 = 0.9;
            const double step = (tau1 - tau0) / static_cast<double>(N_);
            for (int i = 0; i < (N_ + 1); ++i) {
                tau_vec[i] = tau0 + step * static_cast<double>(i);
                const double denom = 1.0 - tau_vec[i];
                r_tau_row[i] = 2.0 / (denom * denom);
            }
        }

        // ---- PRINT: tau and r_tau_row ----
        // print_vec("tau_vec", tau_vec);
        // print_vec("r_tau_row", r_tau_row);

        // -----------------------------
        // 2) Unpack decision vector x
        // -----------------------------
        std::vector<double> px_vector(x, x + (N_ + 1));
        std::vector<double> py_vector(x + 1 * (N_ + 1), x + 2 * (N_ + 1));
        std::vector<double> pz_vector(x + 2 * (N_ + 1), x + 3 * (N_ + 1));
        std::vector<double> psi_vector(x + 3 * (N_ + 1), x + 4 * (N_ + 1));

        std::vector<double> vx_vector(x + 4 * (N_ + 1), x + 5 * (N_ + 1));
        std::vector<double> vy_vector(x + 5 * (N_ + 1), x + 6 * (N_ + 1));
        std::vector<double> vz_vector(x + 6 * (N_ + 1), x + 7 * (N_ + 1));
        std::vector<double> w_vector (x + 7 * (N_ + 1), x + 8 * (N_ + 1));

        std::vector<double> ax_vector(x + 8 * (N_ + 1),  x + 9  * (N_ + 1));
        std::vector<double> ay_vector(x + 9 * (N_ + 1),  x + 10 * (N_ + 1));
        std::vector<double> az_vector(x + 10 * (N_ + 1), x + 11 * (N_ + 1));
        std::vector<double> aw_vector(x + 11 * (N_ + 1), x + 12 * (N_ + 1));

        // ---- PRINT: unpacked decision vectors ----
        // print_vec("px_vector", px_vector);
        // print_vec("py_vector", py_vector);
        // print_vec("pz_vector", pz_vector);
        // print_vec("psi_vector", psi_vector);

        // print_vec("vx_vector", vx_vector);
        // print_vec("vy_vector", vy_vector);
        // print_vec("vz_vector", vz_vector);
        // print_vec("w_vector",  w_vector);

        // print_vec("ax_vector", ax_vector);
        // print_vec("ay_vector", ay_vector);
        // print_vec("az_vector", az_vector);
        // print_vec("aw_vector", aw_vector);

        // -----------------------------
        // 3) Compute tau-derivatives: dyn = Dm^T * X (your current convention)
        // -----------------------------
        std::vector<double> dyn1(N_ + 1), dyn2(N_ + 1), dyn3(N_ + 1), dyn4(N_ + 1);
        std::vector<double> dyn5(N_ + 1), dyn6(N_ + 1), dyn7(N_ + 1), dyn8(N_ + 1);

        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, px_vector.data(), 1, 0.0, dyn1.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, py_vector.data(), 1, 0.0, dyn2.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, pz_vector.data(), 1, 0.0, dyn3.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, psi_vector.data(), 1, 0.0, dyn4.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, vx_vector.data(), 1, 0.0, dyn5.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, vy_vector.data(), 1, 0.0, dyn6.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, vz_vector.data(), 1, 0.0, dyn7.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, w_vector.data(),  1, 0.0, dyn8.data(), 1);

        // ---- PRINT: dyn vectors ----
        // print_vec("dyn1_px_tau",  dyn1);
        // print_vec("dyn2_py_tau",  dyn2);
        // print_vec("dyn3_pz_tau",  dyn3);
        // print_vec("dyn4_psi_tau", dyn4);
        // print_vec("dyn5_vx_tau",  dyn5);
        // print_vec("dyn6_vy_tau",  dyn6);
        // print_vec("dyn7_vz_tau",  dyn7);
        // print_vec("dyn8_w_tau",   dyn8);

        // -----------------------------
        // 9) Apply tau-weighting to RHS
        // -----------------------------
        std::vector<double> px_rhs_scaled(N_ + 1), py_rhs_scaled(N_ + 1), pz_rhs_scaled(N_ + 1), psi_rhs_scaled(N_ + 1);
        std::vector<double> vx_rhs_scaled(N_ + 1), vy_rhs_scaled(N_ + 1), vz_rhs_scaled(N_ + 1), w_rhs_scaled(N_ + 1);

        vdMul(N_ + 1, vx_vector.data(), r_tau_row.data(), px_rhs_scaled.data());
        vdMul(N_ + 1, vy_vector.data(), r_tau_row.data(), py_rhs_scaled.data());
        vdMul(N_ + 1, vz_vector.data(), r_tau_row.data(), pz_rhs_scaled.data());
        vdMul(N_ + 1, w_vector.data(),  r_tau_row.data(), psi_rhs_scaled.data());

        vdMul(N_ + 1, ax_vector.data(), r_tau_row.data(), vx_rhs_scaled.data());
        vdMul(N_ + 1, ay_vector.data(), r_tau_row.data(), vy_rhs_scaled.data());
        vdMul(N_ + 1, az_vector.data(), r_tau_row.data(), vz_rhs_scaled.data());
        vdMul(N_ + 1, aw_vector.data(), r_tau_row.data(), w_rhs_scaled.data());

        // ---- PRINT: RHS scaled ----
        // print_vec("px_rhs_scaled",  px_rhs_scaled);
        // print_vec("py_rhs_scaled",  py_rhs_scaled);
        // print_vec("pz_rhs_scaled",  pz_rhs_scaled);
        // print_vec("psi_rhs_scaled", psi_rhs_scaled);

        // print_vec("vx_rhs_scaled",  vx_rhs_scaled);
        // print_vec("vy_rhs_scaled",  vy_rhs_scaled);
        // print_vec("vz_rhs_scaled",  vz_rhs_scaled);
        // print_vec("w_rhs_scaled",   w_rhs_scaled);

        // -----------------------------
        // 10) Residuals: res = dyn - rhs
        // -----------------------------
        std::vector<double> res_px(N_ + 1), res_py(N_ + 1), res_pz(N_ + 1), res_psi(N_ + 1);
        std::vector<double> res_vx(N_ + 1), res_vy(N_ + 1), res_vz(N_ + 1), res_w(N_ + 1);

        vdSub(N_ + 1, dyn1.data(), px_rhs_scaled.data(),  res_px.data());
        vdSub(N_ + 1, dyn2.data(), py_rhs_scaled.data(),  res_py.data());
        vdSub(N_ + 1, dyn3.data(), pz_rhs_scaled.data(),  res_pz.data());
        vdSub(N_ + 1, dyn4.data(), psi_rhs_scaled.data(), res_psi.data());

        vdSub(N_ + 1, dyn5.data(), vx_rhs_scaled.data(),  res_vx.data());
        vdSub(N_ + 1, dyn6.data(), vy_rhs_scaled.data(),  res_vy.data());
        vdSub(N_ + 1, dyn7.data(), vz_rhs_scaled.data(),  res_vz.data());
        vdSub(N_ + 1, dyn8.data(), w_rhs_scaled.data(),   res_w.data());

        // ---- PRINT: residuals ----
        // print_vec("res_px",  res_px);
        // print_vec("res_py",  res_py);
        // print_vec("res_pz",  res_pz);
        // print_vec("res_psi", res_psi);

        // print_vec("res_vx",  res_vx);
        // print_vec("res_vy",  res_vy);
        // print_vec("res_vz",  res_vz);
        // print_vec("res_w",   res_w);

        // -----------------------------
        // speed and accel constraints
        // -----------------------------
        const double v_max2 = v_max_ * v_max_;
        const double a_max2 = a_max_ * a_max_;
        // print_scalar("v_max2", v_max2);
        // print_scalar("a_max2", a_max2);

        std::vector<double> vx2(N_ + 1), vy2(N_ + 1), vz2(N_ + 1);
        std::vector<double> ax2(N_ + 1), ay2(N_ + 1), az2(N_ + 1);

        vdSqr(N_ + 1, vx_vector.data(), vx2.data());
        vdSqr(N_ + 1, vy_vector.data(), vy2.data());
        vdSqr(N_ + 1, vz_vector.data(), vz2.data());

        vdSqr(N_ + 1, ax_vector.data(), ax2.data());
        vdSqr(N_ + 1, ay_vector.data(), ay2.data());
        vdSqr(N_ + 1, az_vector.data(), az2.data());

        // ---- PRINT: squares ----
        // print_vec("vx2", vx2);
        // print_vec("vy2", vy2);
        // print_vec("vz2", vz2);
        // print_vec("ax2", ax2);
        // print_vec("ay2", ay2);
        // print_vec("az2", az2);

        std::vector<double> v2_sum(N_ + 1), a2_sum(N_ + 1);
        vdAdd(N_ + 1, vx2.data(), vy2.data(), v2_sum.data());
        vdAdd(N_ + 1, v2_sum.data(), vz2.data(), v2_sum.data());

        vdAdd(N_ + 1, ax2.data(), ay2.data(), a2_sum.data());
        vdAdd(N_ + 1, a2_sum.data(), az2.data(), a2_sum.data());

        // ---- PRINT: sums ----
        // print_vec("v2_sum", v2_sum);
        // print_vec("a2_sum", a2_sum);

        std::vector<double> v_max2_vec(N_ + 1, v_max2);
        std::vector<double> a_max2_vec(N_ + 1, a_max2);

        std::vector<double> c_speed(N_ + 1), c_accel(N_ + 1);
        vdSub(N_ + 1, v2_sum.data(), v_max2_vec.data(), c_speed.data());
        vdSub(N_ + 1, a2_sum.data(), a_max2_vec.data(), c_accel.data());

        // ---- PRINT: final inequality constraints ----
        // print_vec("c_speed", c_speed);
        // print_vec("c_accel", c_accel);

        // -----------------------------
        // write g
        // -----------------------------
        for (Index i = 0; i < (N_ + 1); ++i) {
            g[0 * (N_ + 1) + i] = res_px[i];
            g[1 * (N_ + 1) + i] = res_py[i];
            g[2 * (N_ + 1) + i] = res_pz[i];
            g[3 * (N_ + 1) + i] = res_psi[i];
            g[4 * (N_ + 1) + i] = res_vx[i];
            g[5 * (N_ + 1) + i] = res_vy[i];
            g[6 * (N_ + 1) + i] = res_vz[i];
            g[7 * (N_ + 1) + i] = res_w[i];

            g[8 * (N_ + 1) + i] = c_speed[i];
            g[9 * (N_ + 1) + i] = c_accel[i];
        }

        // ---- OPTIONAL: print first few g entries per block ----
        // if (do_print) {
        //     for (int blk = 0; blk < 10; ++blk) {
        //         std::vector<double> tmp(N_ + 1);
        //         for (int i = 0; i < N_ + 1; ++i) tmp[i] = g[blk * (N_ + 1) + i];
        //         print_vec("g_block_" + std::to_string(blk), tmp);
        //     }
        // }

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
        //std::cout << "Finalizing solution" << std::endl;

        //std::cout << "[DEBUG] finalize_solution called. Status = " << status << std::endl;
        //std::cout << "First 10 x: ";
        //for (int i = 0; i < std::min(10, int(n)); ++i) std::cout << x[i] << " ";
        //std::cout << std::endl;


        // Ensure solution_x is 10 * (N + 1)
        solution_x_.resize(12 * (N_ + 1));
        
        // Copy the 10 * (N + 1) elements from the x array
        for (Index i = 0; i < 12 * (N_ + 1); ++i) {
            solution_x_[i] = x[i];
        }
        

        final_obj_value_ = obj_value;         
        bebot_ = Bebot(N_, tf_);
        bebot_.calculate();
        
        // final_time_.resize(1000);
        // for (int i = 0; i < 1000; ++i) {
        //     final_time_[i] = i * tf_ / 999.0;
        // }

        // final_time_.resize(1000);
        // for (int i = 0; i < 1000; ++i) {
        //     final_time_[i] = static_cast<double>(i) / 999.0;  // 0..1
        // }

        const int K = 1000;

        // 1) build BOTH time vectors
        std::vector<double> t_norm(K);
        std::vector<double> t_real(K);

        for (int i = 0; i < K; ++i) {
            const double s = static_cast<double>(i) / static_cast<double>(K - 1); // 0..1
            t_norm[i] = s;          // normalized
            t_real[i] = s * tf_;    // real seconds
        }

        std::vector<double> px_vector(solution_x_.begin(), solution_x_.begin() + (N_ + 1));
        std::vector<double> py_vector(solution_x_.begin() + (N_ + 1), solution_x_.begin() + 2 * (N_ + 1));
        std::vector<double> pz_vector(solution_x_.begin() + 2 * (N_ + 1), solution_x_.begin() + 3 * (N_ + 1));
        std::vector<double> psi_vector(solution_x_.begin() + 3 * (N_ + 1), solution_x_.begin() + 4 * (N_ + 1));
        std::vector<double> vx_vector(solution_x_.begin() + 4 * (N_ + 1), solution_x_.begin() + 5 * (N_ + 1));
        std::vector<double> vy_vector(solution_x_.begin() + 5 * (N_ + 1), solution_x_.begin() + 6 * (N_ + 1));
        std::vector<double> vz_vector(solution_x_.begin() + 6 * (N_ + 1), solution_x_.begin() + 7 * (N_ + 1));
        std::vector<double> w_vector(solution_x_.begin() + 7 * (N_ + 1), solution_x_.begin() + 8 * (N_ + 1));
        std::vector<double> ax_vector(solution_x_.begin() + 8 * (N_ + 1), solution_x_.begin() + 9 * (N_ + 1));
        std::vector<double> ay_vector(solution_x_.begin() + 9 * (N_ + 1), solution_x_.begin() + 10 * (N_ + 1));
        std::vector<double> az_vector(solution_x_.begin() + 10 * (N_ + 1), solution_x_.begin() + 11 * (N_ + 1));
        std::vector<double> aw_vector(solution_x_.begin() + 11 * (N_ + 1), solution_x_.begin() + 12 * (N_ + 1));

        // Helper to print any std::vector<double>
        auto print_vec = [&](const std::string& name, const std::vector<double>& v){
            std::cout << name << " = [";
            for (size_t i = 0; i < v.size(); ++i) {
                std::cout << v[i];
                if (i + 1 < v.size()) std::cout << ", ";
            }
            std::cout << "]\n";
        };

        std::vector<std::vector<double>> px_2d(1, px_vector);
        std::vector<std::vector<double>> py_2d(1, py_vector);
        std::vector<std::vector<double>> pz_2d(1, pz_vector);
        std::vector<std::vector<double>> psi_2d(1, psi_vector);
        std::vector<std::vector<double>> vx_2d(1, vx_vector);
        std::vector<std::vector<double>> vy_2d(1, vy_vector);
        std::vector<std::vector<double>> vz_2d(1, vz_vector);
        std::vector<std::vector<double>> w_2d(1, w_vector);
        std::vector<std::vector<double>> ax_2d(1, ax_vector);
        std::vector<std::vector<double>> ay_2d(1, ay_vector);
        std::vector<std::vector<double>> az_2d(1, az_vector);
        std::vector<std::vector<double>> aw_2d(1, aw_vector);

        // 2) Evaluate Bernstein on normalized domain [0,1]
        auto px_norm = BernsteinPoly(px_2d,     t_norm, 0.0, 1.0);
        auto py_norm = BernsteinPoly(py_2d, t_norm, 0.0, 1.0);
        auto pz_norm = BernsteinPoly(pz_2d,     t_norm, 0.0, 1.0);
        auto psi_norm = BernsteinPoly(psi_2d,     t_norm, 0.0, 1.0);

        auto vx_norm = BernsteinPoly(vx_2d,   t_norm, 0.0, 1.0);
        auto vy_norm = BernsteinPoly(vy_2d,     t_norm, 0.0, 1.0);
        auto vz_norm = BernsteinPoly(vz_2d,     t_norm, 0.0, 1.0);
        auto w_norm = BernsteinPoly(w_2d,     t_norm, 0.0, 1.0);

        auto ax_norm = BernsteinPoly(ax_2d,    t_norm, 0.0, 1.0);
        auto ay_norm = BernsteinPoly(ay_2d,    t_norm, 0.0, 1.0);
        auto az_norm = BernsteinPoly(az_2d,    t_norm, 0.0, 1.0);
        auto aw_norm = BernsteinPoly(aw_2d,    t_norm, 0.0, 1.0);

        // 3) Evaluate Bernstein on real-time domain [0,tf_]
        auto px_real = BernsteinPoly(px_2d,     t_real, 0.0, tf_);
        auto py_real = BernsteinPoly(py_2d, t_real, 0.0, tf_);
        auto pz_real = BernsteinPoly(pz_2d,     t_real, 0.0, tf_);
        auto psi_real = BernsteinPoly(psi_2d,     t_real, 0.0, tf_);

        auto vx_real = BernsteinPoly(vx_2d,   t_real, 0.0, tf_);
        auto vy_real = BernsteinPoly(vy_2d,     t_real, 0.0, tf_);
        auto vz_real = BernsteinPoly(vz_2d,     t_real, 0.0, tf_);
        auto w_real = BernsteinPoly(w_2d,     t_real, 0.0, tf_);

        auto ax_real = BernsteinPoly(ax_2d,    t_real, 0.0, tf_);
        auto ay_real = BernsteinPoly(ay_2d,    t_real, 0.0, tf_);
        auto az_real = BernsteinPoly(az_2d,    t_real, 0.0, tf_);
        auto aw_real = BernsteinPoly(aw_2d,    t_real, 0.0, tf_);

        auto flatten = [](const std::vector<std::vector<double>>& input) {
            std::vector<double> output;
            for (const auto& row : input) {
                output.insert(output.end(), row.begin(), row.end());
            }
            return output;
        };
        
        // ---- real time outputs (0..tf_) ----
        // States: p and psi
        writeToCSV(t_real, flatten(px_real),  "px_real.csv");
        writeToCSV(bebot_.getNodes(), px_vector, "px_controlpoints_real.csv");

        writeToCSV(t_real, flatten(py_real),  "py_real.csv");
        writeToCSV(bebot_.getNodes(), py_vector, "py_controlpoints_real.csv");

        writeToCSV(t_real, flatten(pz_real),  "pz_real.csv");
        writeToCSV(bebot_.getNodes(), pz_vector, "pz_controlpoints_real.csv");

        writeToCSV(t_real, flatten(psi_real), "psi_real.csv");
        writeToCSV(bebot_.getNodes(), psi_vector, "psi_controlpoints_real.csv");

        // Velocities: v and w
        writeToCSV(t_real, flatten(vx_real),  "vx_real.csv");
        writeToCSV(bebot_.getNodes(), vx_vector, "vx_controlpoints_real.csv");

        writeToCSV(t_real, flatten(vy_real),  "vy_real.csv");
        writeToCSV(bebot_.getNodes(), vy_vector, "vy_controlpoints_real.csv");

        writeToCSV(t_real, flatten(vz_real),  "vz_real.csv");
        writeToCSV(bebot_.getNodes(), vz_vector, "vz_controlpoints_real.csv");

        writeToCSV(t_real, flatten(w_real),   "w_real.csv");
        writeToCSV(bebot_.getNodes(), w_vector, "w_controlpoints_real.csv");

        // Controls: a and aw
        writeToCSV(t_real, flatten(ax_real),  "ax_real.csv");
        writeToCSV(bebot_.getNodes(), ax_vector, "ax_controlpoints_real.csv");

        writeToCSV(t_real, flatten(ay_real),  "ay_real.csv");
        writeToCSV(bebot_.getNodes(), ay_vector, "ay_controlpoints_real.csv");

        writeToCSV(t_real, flatten(az_real),  "az_real.csv");
        writeToCSV(bebot_.getNodes(), az_vector, "az_controlpoints_real.csv");

        writeToCSV(t_real, flatten(aw_real),  "aw_real.csv");
        writeToCSV(bebot_.getNodes(), aw_vector, "aw_controlpoints_real.csv");

    }

    const std::vector<Number>& get_solution_x() const { return solution_x_; }
    Number get_final_obj_value() const { return final_obj_value_; }

private:
    int N_;
    double tf_;
    double px_max_;
    double px_min_;
    double py_max_;
    double py_min_;
    double pz_max_;
    double pz_min_;
    double psi_max_;
    double psi_min_;
    double v_max_;
    double w_max_;
    double a_max_;
    double aw_max_;
    double px_cur_;
    double py_cur_;
    double pz_cur_;
    double psi_cur_;
    double vx_cur_;
    double vy_cur_;
    double vz_cur_;
    double w_cur_;
    double pxf_;
    double pyf_;
    double pzf_;
    double psif_;
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

extern "C" {
    PointSetProblem* create_point_set_problem(int N, double tf, 
        double px_max, double px_min, double py_max, double py_min, double pz_max, double pz_min, double psi_max, double psi_min,   
        double v_max, double w_max, 
        double a_max, double aw_max, 
        double px_cur, double py_cur, double pz_cur, double psi_cur,
        double vx_cur, double vy_cur, double vz_cur, double w_cur,
        double pxf, double pyf, double pzf, double psif) {
        
        return new PointSetProblem(N, tf, 
            px_max, px_min, py_max, py_min, pz_max, pz_min, psi_max, psi_min, 
            v_max, w_max, 
            a_max, aw_max,
            px_cur, py_cur, pz_cur, psi_cur,
            vx_cur, vy_cur, vz_cur, w_cur,
            pxf, pyf, pzf, psif);
    }

    void solve_point_set_problem(PointSetProblem* problem) {
        SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
        app->Options()->SetStringValue("linear_solver", "ma57");
        app->Options()->SetStringValue("mu_strategy", "adaptive");
        app->Options()->SetStringValue("gradient_approximation", "finite-difference-values");
        app->Options()->SetStringValue("jacobian_approximation", "finite-difference-values");
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");
        app->Options()->SetIntegerValue("max_iter", 1000);
        app->Options()->SetNumericValue("tol",             1e-6);   // OptimalityTolerance = 1e-3    1
        app->Options()->SetNumericValue("constr_viol_tol", 1e-6);
        app->Options()->SetNumericValue("obj_scaling_factor", 1e-6);
    

        //app->Options()->SetNumericValue("constr_viol_tol", 1e-6);
        app->Options()->SetIntegerValue("print_level", 5); 
        //app->Options()->SetStringValue("nlp_scaling_method", "gradient-based");


        app->RethrowNonIpoptException(true);
        ApplicationReturnStatus status = app->Initialize();
        if (status != Solve_Succeeded) {
            std::cerr << "IPOPT initialization failed!" << std::endl;
            return;
        }
        status = app->OptimizeTNLP(problem);
        if (status == Solve_Succeeded || status == Solved_To_Acceptable_Level) {
            std::cout << "Optimization succeeded!" << std::endl;
            
            const auto& solution_x = problem->get_solution_x();
            std::cout << "Optimal Solution (x): ";
            for (Index i = 0; i < solution_x.size(); i++) {
                std::cout << solution_x[i] << " ";
            }
            std::cout << std::endl;
        } else {
            std::cerr << "Optimization failed with status " << status << std::endl;
        }
    }

    void get_solution(PointSetProblem* problem, double* solution, int n) {
        const std::vector<double>& sol = problem->get_solution_x();

        std::cout << "[DEBUG] get_solution called. Vector size: " << sol.size() << std::endl;
        // Print all values in one line, comma-separated
        //std::cout << "[DEBUG] Full solution vector: [";
        //for (size_t i = 0; i < sol.size(); ++i) {
        //    std::cout << sol[i];
        //    if (i + 1 < sol.size()) std::cout << ", ";
        //}
        //std::cout << "]" << std::endl;

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

// g++ -shared -fPIC -o libbebot_mpc_drone_single.so ~/dev/optimization/BeBOT_cpp_v2/examples/bebot/ma_57/example_mpc_drone_threed/libbebot_mpc_drone_single.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bebot.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinpoly.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteindifferentialmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinmatrix_a2b.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/degelevmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/nchoosek_mod.cpp -I~/dev/optimization/BeBOT_cpp_v2/include -I./Ipopt/src/ -L./Ipopt/src/.libs -lipopt -L/opt/intel/oneapi/mkl/latest/lib/intel64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -ldl -lm -lpthread -lstdc++
