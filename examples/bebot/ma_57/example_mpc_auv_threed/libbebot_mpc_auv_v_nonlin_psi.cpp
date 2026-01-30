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
#include "state_space_matrices.h"

#include <array>
#include <string>
#include <limits>
#include <algorithm>

using namespace Ipopt;

class PointSetProblem : public Ipopt::TNLP {
public:
    PointSetProblem(int N, double tf, double delta_v_max, double delta_v_min,
        double delta_s_max, double delta_s_min, double delta_m_max, double delta_m_min, double delta_h_max, double delta_h_min,
        double zmax, double zmin, double wmax, double wmin, double thetamax, double thetamin, double qmax, double qmin,
        double psimax, double psimin, double rmax, double rmin,
        double z0, double w0, double theta0, double q0,
        double psi0,  double r0,
        double delta_v0, double delta_s0, double delta_m0, double delta_h0,
        double zf, double thetaf,
        double psif,
        double a11, double a12, double a13, double a14,
        double a21, double a22, double a23, double a24,
        double a31, double a32, double a33, double a34,
        double a41, double a42, double a43, double a44,
        double b11, double b12, double b13,
        double b21, double b22, double b23,
        double b31, double b32, double b33,
        double b41, double b42, double b43,
        double c11, double c12,
        double c21, double c22,
        double d11,
        double d21,
        double t0, double tend,
        double dv_prev, double ds_prev, double dm_prev, double dh_prev,
        double dv_dot_max, double dm_dot_max, double ds_dot_max, double dh_dot_max)
        : N_(N), tf_(tf), delta_v_max_(delta_v_max), delta_v_min_(delta_v_min), delta_m_max_(delta_m_max), delta_m_min_(delta_m_min),
        delta_s_max_(delta_s_max), delta_s_min_(delta_s_min), delta_h_max_(delta_h_max), delta_h_min_(delta_h_min),
        zmax_(zmax), zmin_(zmin), wmax_(wmax), wmin_(wmin), thetamax_(thetamax), thetamin_(thetamin), qmax_(qmax), qmin_(qmin),
        psimax_(psimax), psimin_(psimin), rmax_(rmax), rmin_(rmin),
        z0_(z0), w0_(w0), theta0_(theta0), q0_(q0),
        psi0_(psi0), r0_(r0),
        delta_v0_(delta_v0), delta_s0_(delta_s0), delta_m0_(delta_m0), delta_h0_(delta_h0),
        zf_(zf), thetaf_(thetaf),
        psif_(psif), bebot_(N, tf_),
        t0_(t0), tend_(tend),
        dv_prev_(dv_prev), ds_prev_(ds_prev), dm_prev_(dm_prev), dh_prev_(dh_prev),
        dv_dot_max_(dv_dot_max), dm_dot_max_(dm_dot_max), ds_dot_max_(ds_dot_max), dh_dot_max_(dh_dot_max) {

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

        C_ = {{
            {c11, c12},
            {c21, c22}
        }};

        D_ = {{
            {d11},
            {d21}
        }};

        bebot_.calculate();

        // ====== CACHE “COMPUTE-ONCE” STUFF HERE (no math changes) ======
        build_caches_();
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
        n = 10 * (N_ + 1); // z, theta, w, q, psi, r, dv, dm, ds, dh
        m = 10 * (N_ + 1); // 6*(N+1) dynamics + 4*(N+1) control derivatives
        nnz_jac_g = n * m;
        nnz_h_lag = 0;
        index_style = TNLP::C_STYLE;
        return true;
    }

    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u) {
        std::vector<double> x_lower(n, -std::numeric_limits<double>::infinity());
        std::vector<double> x_upper(n, std::numeric_limits<double>::infinity());

        // z
        for (int i = 1; i < N_ + 1; ++i) { x_lower[i] = zmin_; x_upper[i] = zmax_; }
        x_lower[0] = x_upper[0] = z0_;

        // theta
        for (int i = (N_ + 1) + 1; i < 2 * (N_ + 1); ++i) { x_lower[i] = thetamin_; x_upper[i] = thetamax_; }
        x_lower[(N_ + 1)] = x_upper[(N_ + 1)] = theta0_;

        // w
        for (int i = 2 * (N_ + 1) + 1; i < 3 * (N_ + 1); ++i) { x_lower[i] = wmin_; x_upper[i] = wmax_; }
        x_lower[2 * (N_ + 1)] = x_upper[2 * (N_ + 1)] = w0_;

        // q
        for (int i = 3 * (N_ + 1) + 1; i < 4 * (N_ + 1); ++i) { x_lower[i] = qmin_; x_upper[i] = qmax_; }
        x_lower[3 * (N_ + 1)] = x_upper[3 * (N_ + 1)] = q0_;

        // psi
        for (int i = 4 * (N_ + 1) + 1; i < 5 * (N_ + 1); ++i) { x_lower[i] = psimin_; x_upper[i] = psimax_; }
        x_lower[4 * (N_ + 1)] = x_upper[4 * (N_ + 1)] = psi0_;

        // r
        for (int i = 5 * (N_ + 1) + 1; i < 6 * (N_ + 1); ++i) { x_lower[i] = rmin_; x_upper[i] = rmax_; }
        x_lower[5 * (N_ + 1)] = x_upper[5 * (N_ + 1)] = r0_;

        // dv
        for (int i = 6 * (N_ + 1); i < 7 * (N_ + 1); ++i) { x_lower[i] = delta_v_min_; x_upper[i] = delta_v_max_; }
        x_lower[6 * (N_ + 1)] = x_upper[6 * (N_ + 1)] = dv_prev_;

        // dm
        for (int i = 7 * (N_ + 1); i < 8 * (N_ + 1); ++i) { x_lower[i] = delta_m_min_; x_upper[i] = delta_m_max_; }
        x_lower[7 * (N_ + 1)] = x_upper[7 * (N_ + 1)] = dm_prev_;

        // ds
        for (int i = 8 * (N_ + 1); i < 9 * (N_ + 1); ++i) { x_lower[i] = delta_s_min_; x_upper[i] = delta_s_max_; }
        x_lower[8 * (N_ + 1)] = x_upper[8 * (N_ + 1)] = ds_prev_;

        // dh
        for (int i = 9 * (N_ + 1); i < 10 * (N_ + 1); ++i) { x_lower[i] = delta_h_min_; x_upper[i] = delta_h_max_; }
        x_lower[9 * (N_ + 1)] = x_upper[9 * (N_ + 1)] = dh_prev_;

        std::copy(x_lower.begin(), x_lower.end(), x_l);
        std::copy(x_upper.begin(), x_upper.end(), x_u);

        // dynamics equalities: 6*(N+1) constraints -> g == 0
        for (int i = 0; i < 6 * (N_ + 1); ++i) { g_l[i] = 0.0; g_u[i] = 0.0; }

        // dv_dt bounds
        for (int i = 6 * (N_ + 1); i < 7 * (N_ + 1); ++i) { g_l[i] = -dv_dot_max_; g_u[i] = dv_dot_max_; }

        // dm_dt bounds
        for (int i = 7 * (N_ + 1); i < 8 * (N_ + 1); ++i) { g_l[i] = -dm_dot_max_; g_u[i] = dm_dot_max_; }

        // ds_dt bounds
        for (int i = 8 * (N_ + 1); i < 9 * (N_ + 1); ++i) { g_l[i] = -ds_dot_max_; g_u[i] = ds_dot_max_; }

        // dh_dt bounds
        for (int i = 9 * (N_ + 1); i < 10 * (N_ + 1); ++i) { g_l[i] = -dh_dot_max_; g_u[i] = dh_dot_max_; }

        return true;
    }

    virtual bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda) {
        for (Index i = 0; i < n; ++i) { x[i] = 1.0; }
        return true;
    }

    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value) override
    {
        // ===== weights (unchanged) =====
        const double wz     = 100;
        const double wtheta = 1;
        const double wpsi   = 10;

        const double wdv = 100;
        const double wdm = 0.01;
        const double wds = 0.01;
        const double wdh = 0.001;

        const int Np1 = N_ + 1;

        // ===== workspace: allocate once per thread, reuse every call =====
        struct WS {
            int Np1 = -1;
            double zf = std::numeric_limits<double>::quiet_NaN();
            double thetaf = std::numeric_limits<double>::quiet_NaN();
            double psif = std::numeric_limits<double>::quiet_NaN();

            std::vector<double> zf_vec, thetaf_vec, psif_vec;
            std::vector<double> r_tau, ones;

            std::vector<double> z_diff, theta_diff, psi_diff;
            std::vector<double> z_diff_sqr, theta_diff_sqr, psi_diff_sqr;
            std::vector<double> z_weighted, theta_weighted, psi_weighted;

            std::vector<double> dv_sqr, dm_sqr, ds_sqr, dh_sqr;

            void ensure(int Np1_in, double zf_in, double thetaf_in, double psif_in)
            {
                if (Np1 != Np1_in) {
                    Np1 = Np1_in;

                    zf_vec.resize(Np1);
                    thetaf_vec.resize(Np1);
                    psif_vec.resize(Np1);

                    r_tau.resize(Np1);
                    ones.assign(Np1, 1.0);

                    z_diff.resize(Np1);
                    theta_diff.resize(Np1);
                    psi_diff.resize(Np1);

                    z_diff_sqr.resize(Np1);
                    theta_diff_sqr.resize(Np1);
                    psi_diff_sqr.resize(Np1);

                    z_weighted.resize(Np1);
                    theta_weighted.resize(Np1);
                    psi_weighted.resize(Np1);

                    dv_sqr.resize(Np1);
                    dm_sqr.resize(Np1);
                    ds_sqr.resize(Np1);
                    dh_sqr.resize(Np1);

                    // r_tau = 2 ./ (1 - tau).^2, tau = linspace(-1, 0.999, N+1)
                    if (Np1 == 1) {
                        const double tau0 = -1.0;
                        const double denom = (1.0 - tau0);
                        r_tau[0] = 2.0 / (denom * denom);
                    } else {
                        const double tau0 = -1.0;
                        const double tau1 = 0.999;
                        const double step = (tau1 - tau0) / static_cast<double>(Np1 - 1);
                        for (int i = 0; i < Np1; ++i) {
                            const double tau = tau0 + step * static_cast<double>(i);
                            const double denom = (1.0 - tau);
                            r_tau[i] = 2.0 / (denom * denom);
                        }
                    }
                }

                // (re)fill target vectors only if targets changed
                if (zf != zf_in) {
                    zf = zf_in;
                    std::fill(zf_vec.begin(), zf_vec.end(), zf);
                }
                if (thetaf != thetaf_in) {
                    thetaf = thetaf_in;
                    std::fill(thetaf_vec.begin(), thetaf_vec.end(), thetaf);
                }
                if (psif != psif_in) {
                    psif = psif_in;
                    std::fill(psif_vec.begin(), psif_vec.end(), psif);
                }
            }
        };
        static thread_local WS ws;
        ws.ensure(Np1, zf_, thetaf_, psif_);

        // x layout: [z, theta, w, q, psi, r, dv, dm, ds, dh], each length (N+1)
        const double* z_ptr     = reinterpret_cast<const double*>(x) + 0 * Np1;
        const double* theta_ptr = reinterpret_cast<const double*>(x) + 1 * Np1;
        const double* psi_ptr   = reinterpret_cast<const double*>(x) + 4 * Np1;

        const double* dv_ptr = reinterpret_cast<const double*>(x) + 6 * Np1;
        const double* dm_ptr = reinterpret_cast<const double*>(x) + 7 * Np1;
        const double* ds_ptr = reinterpret_cast<const double*>(x) + 8 * Np1;
        const double* dh_ptr = reinterpret_cast<const double*>(x) + 9 * Np1;

        // --- errors ---
        vdSub(Np1, z_ptr,     ws.zf_vec.data(),     ws.z_diff.data());
        vdSub(Np1, theta_ptr, ws.thetaf_vec.data(), ws.theta_diff.data());
        vdSub(Np1, psi_ptr,   ws.psif_vec.data(),   ws.psi_diff.data());

        // --- squared errors ---
        vdSqr(Np1, ws.z_diff.data(),     ws.z_diff_sqr.data());
        vdSqr(Np1, ws.theta_diff.data(), ws.theta_diff_sqr.data());
        vdSqr(Np1, ws.psi_diff.data(),   ws.psi_diff_sqr.data());

        // --- weighted squared errors ---
        vdMul(Np1, ws.r_tau.data(), ws.z_diff_sqr.data(),     ws.z_weighted.data());
        vdMul(Np1, ws.r_tau.data(), ws.theta_diff_sqr.data(), ws.theta_weighted.data());
        vdMul(Np1, ws.r_tau.data(), ws.psi_diff_sqr.data(),   ws.psi_weighted.data());

        // sums: arrays are nonnegative (r_tau>0, squares>=0) so dasum == sum
        const double sum_z_state     = cblas_dasum(Np1, ws.z_weighted.data(),     1);
        const double sum_theta_state = cblas_dasum(Np1, ws.theta_weighted.data(), 1);
        const double sum_psi_state   = cblas_dasum(Np1, ws.psi_weighted.data(),   1);

        const double state_term =
            wz     * sum_z_state +
            wtheta * sum_theta_state +
            wpsi   * sum_psi_state;

        // --- control term ---
        vdSqr(Np1, dv_ptr, ws.dv_sqr.data());
        vdSqr(Np1, dm_ptr, ws.dm_sqr.data());
        vdSqr(Np1, ds_ptr, ws.ds_sqr.data());
        vdSqr(Np1, dh_ptr, ws.dh_sqr.data());

        const double sum_dv = cblas_dasum(Np1, ws.dv_sqr.data(), 1);
        const double sum_dm = cblas_dasum(Np1, ws.dm_sqr.data(), 1);
        const double sum_ds = cblas_dasum(Np1, ws.ds_sqr.data(), 1);
        const double sum_dh = cblas_dasum(Np1, ws.dh_sqr.data(), 1);

        const double control_term =
            wdv * sum_dv +
            wdm * sum_dm +
            wds * sum_ds +
            wdh * sum_dh;

        obj_value = state_term + control_term;
        return true;
    }


    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {

        // ====== cached Dm, tau_vec, r_tau_row (same as before) ======
        const std::vector<double>& Dm       = Dm_cached_;
        const std::vector<double>& tau_vec  = tau_vec_cached_;
        const std::vector<double>& r_tau_row = r_tau_row_cached_;

        // -----------------------------
        // 2) Unpack decision vector x into state/control knot vectors
        // -----------------------------
        std::vector<double> z_vector(x, x + (N_ + 1));
        std::vector<double> theta_vector(x + 1 * (N_ + 1), x + 2 * (N_ + 1));
        std::vector<double> w_vector(x + 2 * (N_ + 1), x + 3 * (N_ + 1));
        std::vector<double> q_vector(x + 3 * (N_ + 1), x + 4 * (N_ + 1));
        std::vector<double> psi_vector(x + 4 * (N_ + 1), x + 5 * (N_ + 1));
        std::vector<double> r_vector(x + 5 * (N_ + 1), x + 6 * (N_ + 1));

        std::vector<double> delta_v_vector(x + 6 * (N_ + 1), x + 7 * (N_ + 1));
        std::vector<double> delta_m_vector(x + 7 * (N_ + 1), x + 8 * (N_ + 1));
        std::vector<double> delta_s_vector(x + 8 * (N_ + 1), x + 9 * (N_ + 1));
        std::vector<double> delta_h_vector(x + 9 * (N_ + 1), x + 10 * (N_ + 1));

        // -----------------------------
        // 3) Compute tau-derivatives of states via BeBOT: X_tau = X * Dm
        // -----------------------------
        std::vector<double> dyn1(N_ + 1);
        std::vector<double> dyn2(N_ + 1);
        std::vector<double> dyn3(N_ + 1);
        std::vector<double> dyn4(N_ + 1);
        std::vector<double> dyn5(N_ + 1);
        std::vector<double> dyn6(N_ + 1);

        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, z_vector.data(), 1, 0.0, dyn1.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, theta_vector.data(), 1, 0.0, dyn2.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, w_vector.data(), 1, 0.0, dyn3.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, q_vector.data(), 1, 0.0, dyn4.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, psi_vector.data(), 1, 0.0, dyn5.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, r_vector.data(), 1, 0.0, dyn6.data(), 1);

        // -----------------------------
        // 4) Build depth subsystem matrices Xd (4x(N+1)) and Ud (3x(N+1))
        // -----------------------------
        std::vector<double> X1_matrix_flat(4 * (N_ + 1));
        std::vector<double> U_matrix_flat(3 * (N_ + 1));

        for (Index i = 0; i < (N_ + 1); ++i) {
            X1_matrix_flat[i]                 = z_vector[i];
            X1_matrix_flat[(N_ + 1) + i]      = theta_vector[i];
            X1_matrix_flat[2 * (N_ + 1) + i]  = w_vector[i];
            X1_matrix_flat[3 * (N_ + 1) + i]  = q_vector[i];
        }

        for (Index i = 0; i < (N_ + 1); ++i) {
            U_matrix_flat[i]                 = delta_v_vector[i];
            U_matrix_flat[(N_ + 1) + i]      = delta_m_vector[i];
            U_matrix_flat[2 * (N_ + 1) + i]  = delta_s_vector[i];
        }

        // -----------------------------
        // 5) Compute depth RHS: (A*Xd + B*Ud)
        // -----------------------------
        std::vector<double> X2_matrix_flat(4 * (N_ + 1), 0.0);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4, (N_ + 1), 4, 1.0,
                    &A_[0][0], 4, X1_matrix_flat.data(), (N_ + 1),
                    0.0, X2_matrix_flat.data(), (N_ + 1));

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4, (N_ + 1), 3, 1.0,
                    &B_[0][0], 3, U_matrix_flat.data(), (N_ + 1),
                    1.0, X2_matrix_flat.data(), (N_ + 1));

        // -----------------------------
        // 6) Build horizontal subsystem matrices Xh (2x(N+1)) and Uh (1x(N+1))
        // -----------------------------
        std::vector<double> X1_matrix_flat_cd(2 * (N_ + 1));
        std::vector<double> U_matrix_flat_cd(1 * (N_ + 1));

        for (Index i = 0; i < (N_ + 1); ++i) {
            X1_matrix_flat_cd[i]            = psi_vector[i];
            X1_matrix_flat_cd[(N_ + 1) + i] = r_vector[i];
        }
        for (Index i = 0; i < (N_ + 1); ++i) {
            U_matrix_flat_cd[i] = delta_h_vector[i];
        }

        // -----------------------------
        // 7) Compute horizontal RHS: (C*Xh + D*Uh)
        // -----------------------------
        std::vector<double> X2_matrix_flat_cd(2 * (N_ + 1), 0.0);

        // (your original file had C*Xh commented out; left as-is)
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, (N_ + 1), 2, 1.0, &C_[0][0], 2, X1_matrix_flat_cd.data(), (N_ + 1),0.0, X2_matrix_flat_cd.data(), (N_ + 1));
        // pm("X2_cd after C*Xh", X2_matrix_flat_cd.data(), 2, (N_+1), /*rowMajor=*/true);
        // X2_cd += D*Uh
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, (N_ + 1), 1, 1.0, &D_[0][0], 1, U_matrix_flat_cd.data(), (N_ + 1),1.0, X2_matrix_flat_cd.data(), (N_ + 1));

        // -----------------------------
        // 8) Unpack RHS rows into per-state vectors
        // -----------------------------
        std::vector<double> z_rhs(N_ + 1);
        std::vector<double> theta_rhs(N_ + 1);
        std::vector<double> w_rhs(N_ + 1);
        std::vector<double> q_rhs(N_ + 1);

        std::vector<double> psi_rhs(N_ + 1);
        std::vector<double> r_rhs(N_ + 1);

        for (Index i = 0; i < (N_ + 1); ++i) {
            z_rhs[i]     = X2_matrix_flat[i];
            theta_rhs[i] = X2_matrix_flat[(N_ + 1) + i];
            w_rhs[i]     = X2_matrix_flat[2 * (N_ + 1) + i];
            q_rhs[i]     = X2_matrix_flat[3 * (N_ + 1) + i];
        }
        for (Index i = 0; i < (N_ + 1); ++i) {
            psi_rhs[i] = X2_matrix_flat_cd[i];
            r_rhs[i]   = X2_matrix_flat_cd[(N_ + 1) + i];
        }

        // -----------------------------
        // 9) Apply tau-weighting to RHS
        // -----------------------------
        std::vector<double> z_rhs_scaled(N_ + 1);
        std::vector<double> theta_rhs_scaled(N_ + 1);
        std::vector<double> w_rhs_scaled(N_ + 1);
        std::vector<double> q_rhs_scaled(N_ + 1);

        std::vector<double> psi_rhs_scaled(N_ + 1);
        std::vector<double> r_rhs_scaled(N_ + 1);

        vdMul(N_ + 1, z_rhs.data(),     r_tau_row.data(), z_rhs_scaled.data());
        vdMul(N_ + 1, theta_rhs.data(), r_tau_row.data(), theta_rhs_scaled.data());
        vdMul(N_ + 1, w_rhs.data(),     r_tau_row.data(), w_rhs_scaled.data());
        vdMul(N_ + 1, q_rhs.data(),     r_tau_row.data(), q_rhs_scaled.data());

        vdMul(N_ + 1, psi_rhs.data(),   r_tau_row.data(), psi_rhs_scaled.data());
        vdMul(N_ + 1, r_rhs.data(),     r_tau_row.data(), r_rhs_scaled.data());

        // -----------------------------
        // 10) Dynamics equality residuals
        // -----------------------------
        std::vector<double> res_z(N_ + 1);
        std::vector<double> res_theta(N_ + 1);
        std::vector<double> res_w(N_ + 1);
        std::vector<double> res_q(N_ + 1);
        std::vector<double> res_psi(N_ + 1);
        std::vector<double> res_r(N_ + 1);

        vdSub(N_ + 1, dyn1.data(), z_rhs_scaled.data(),     res_z.data());
        vdSub(N_ + 1, dyn2.data(), theta_rhs_scaled.data(), res_theta.data());
        vdSub(N_ + 1, dyn3.data(), w_rhs_scaled.data(),     res_w.data());
        vdSub(N_ + 1, dyn4.data(), q_rhs_scaled.data(),     res_q.data());
        vdSub(N_ + 1, dyn5.data(), psi_rhs_scaled.data(),   res_psi.data());
        vdSub(N_ + 1, dyn6.data(), r_rhs_scaled.data(),     res_r.data());

        // -----------------------------
        // 11) Control derivatives
        // -----------------------------
        std::vector<double> dyn7(N_ + 1);
        std::vector<double> dyn8(N_ + 1);
        std::vector<double> dyn9(N_ + 1);
        std::vector<double> dyn10(N_ + 1);

        std::vector<double> dv_tau(N_ + 1);
        std::vector<double> dm_tau(N_ + 1);
        std::vector<double> ds_tau(N_ + 1);
        std::vector<double> dh_tau(N_ + 1);

        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, delta_v_vector.data(), 1, 0.0, dv_tau.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, delta_m_vector.data(), 1, 0.0, dm_tau.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, delta_s_vector.data(), 1, 0.0, ds_tau.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, delta_h_vector.data(), 1, 0.0, dh_tau.data(), 1);

        vdDiv(N_ + 1, dv_tau.data(), r_tau_row.data(), dyn7.data());
        vdDiv(N_ + 1, dm_tau.data(), r_tau_row.data(), dyn8.data());
        vdDiv(N_ + 1, ds_tau.data(), r_tau_row.data(), dyn9.data());
        vdDiv(N_ + 1, dh_tau.data(), r_tau_row.data(), dyn10.data());

        // -----------------------------
        // 12) Pack constraints into g
        // -----------------------------
        for (Index i = 0; i < (N_ + 1); ++i) {
            g[0 * (N_ + 1) + i] = res_z[i];
            g[1 * (N_ + 1) + i] = res_theta[i];
            g[2 * (N_ + 1) + i] = res_w[i];
            g[3 * (N_ + 1) + i] = res_q[i];
            g[4 * (N_ + 1) + i] = res_psi[i];
            g[5 * (N_ + 1) + i] = res_r[i];

            g[6 * (N_ + 1) + i] = dyn7[i];
            g[7 * (N_ + 1) + i] = dyn8[i];
            g[8 * (N_ + 1) + i] = dyn9[i];
            g[9 * (N_ + 1) + i] = dyn10[i];
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
        solution_x_.resize(10 * (N_ + 1));
        for (Index i = 0; i < 10 * (N_ + 1); ++i) {
            solution_x_[i] = x[i];
        }

        final_obj_value_ = obj_value;
        bebot_ = Bebot(N_, tf_);
        bebot_.calculate();

        const int K = 1000;
        std::vector<double> t_norm(K);
        std::vector<double> t_real(K);

        for (int i = 0; i < K; ++i) {
            const double s = static_cast<double>(i) / static_cast<double>(K - 1);
            t_norm[i] = s;
            t_real[i] = s * tf_;
        }

        std::vector<double> z_vector(solution_x_.begin(), solution_x_.begin() + (N_ + 1));
        std::vector<double> theta_vector(solution_x_.begin() + (N_ + 1), solution_x_.begin() + 2 * (N_ + 1));
        std::vector<double> w_vector(solution_x_.begin() + 2 * (N_ + 1), solution_x_.begin() + 3 * (N_ + 1));
        std::vector<double> q_vector(solution_x_.begin() + 3 * (N_ + 1), solution_x_.begin() + 4 * (N_ + 1));
        std::vector<double> psi_vector(solution_x_.begin() + 4 * (N_ + 1), solution_x_.begin() + 5 * (N_ + 1));
        std::vector<double> r_vector(solution_x_.begin() + 5 * (N_ + 1), solution_x_.begin() + 6 * (N_ + 1));
        std::vector<double> delta_v_vector(solution_x_.begin() + 6 * (N_ + 1), solution_x_.begin() + 7 * (N_ + 1));
        std::vector<double> delta_m_vector(solution_x_.begin() + 7 * (N_ + 1), solution_x_.begin() + 8 * (N_ + 1));
        std::vector<double> delta_s_vector(solution_x_.begin() + 8 * (N_ + 1), solution_x_.begin() + 9 * (N_ + 1));
        std::vector<double> delta_h_vector(solution_x_.begin() + 9 * (N_ + 1), solution_x_.begin() + 10 * (N_ + 1));

        auto flatten = [](const std::vector<std::vector<double>>& input) {
            std::vector<double> output;
            for (const auto& row : input) {
                output.insert(output.end(), row.begin(), row.end());
            }
            return output;
        };

        std::vector<std::vector<double>> z_2d(1, z_vector);
        std::vector<std::vector<double>> theta_2d(1, theta_vector);
        std::vector<std::vector<double>> w_2d(1, w_vector);
        std::vector<std::vector<double>> q_2d(1, q_vector);
        std::vector<std::vector<double>> psi_2d(1, psi_vector);
        std::vector<std::vector<double>> r_2d(1, r_vector);
        std::vector<std::vector<double>> delta_v_2d(1, delta_v_vector);
        std::vector<std::vector<double>> delta_m_2d(1, delta_m_vector);
        std::vector<std::vector<double>> delta_s_2d(1, delta_s_vector);
        std::vector<std::vector<double>> delta_h_2d(1, delta_h_vector);

        auto z_real     = BernsteinPoly(z_2d,     t_real, 0.0, tf_);
        auto theta_real = BernsteinPoly(theta_2d, t_real, 0.0, tf_);
        auto w_real     = BernsteinPoly(w_2d,     t_real, 0.0, tf_);
        auto q_real     = BernsteinPoly(q_2d,     t_real, 0.0, tf_);
        auto psi_real   = BernsteinPoly(psi_2d,   t_real, 0.0, tf_);
        auto r_real     = BernsteinPoly(r_2d,     t_real, 0.0, tf_);
        auto dv_real    = BernsteinPoly(delta_v_2d, t_real, 0.0, tf_);
        auto dm_real    = BernsteinPoly(delta_m_2d, t_real, 0.0, tf_);
        auto ds_real    = BernsteinPoly(delta_s_2d, t_real, 0.0, tf_);
        auto dh_real    = BernsteinPoly(delta_h_2d, t_real, 0.0, tf_);

        writeToCSV(t_real, flatten(z_real),     "z_real.csv");
        writeToCSV(bebot_.getNodes(), z_vector, "z_controlpoints_real.csv");

        writeToCSV(t_real, flatten(w_real),     "w_real.csv");
        writeToCSV(bebot_.getNodes(), w_vector, "w_controlpoints_real.csv");

        writeToCSV(t_real, flatten(theta_real),     "theta_real.csv");
        writeToCSV(bebot_.getNodes(), theta_vector, "theta_controlpoints_real.csv");

        writeToCSV(t_real, flatten(q_real),     "q_real.csv");
        writeToCSV(bebot_.getNodes(), q_vector, "q_controlpoints_real.csv");

        writeToCSV(t_real, flatten(psi_real),     "psi_real.csv");
        writeToCSV(bebot_.getNodes(), psi_vector, "psi_controlpoints_real.csv");

        writeToCSV(t_real, flatten(r_real),     "r_real.csv");
        writeToCSV(bebot_.getNodes(), r_vector, "r_controlpoints_real.csv");

        writeToCSV(t_real, flatten(dv_real),          "delta_v_real.csv");
        writeToCSV(bebot_.getNodes(), delta_v_vector, "delta_v_controlpoints_real.csv");

        writeToCSV(t_real, flatten(dm_real),          "delta_m_real.csv");
        writeToCSV(bebot_.getNodes(), delta_m_vector, "delta_m_controlpoints_real.csv");

        writeToCSV(t_real, flatten(ds_real),          "delta_s_real.csv");
        writeToCSV(bebot_.getNodes(), delta_s_vector, "delta_s_controlpoints_real.csv");

        writeToCSV(t_real, flatten(dh_real),          "delta_h_real.csv");
        writeToCSV(bebot_.getNodes(), delta_h_vector, "delta_h_controlpoints_real.csv");
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
    double delta_h_max_;
    double delta_h_min_;
    double zmax_;
    double zmin_;
    double wmax_;
    double wmin_;
    double thetamax_;
    double thetamin_;
    double qmax_;
    double qmin_;
    double psimax_;
    double psimin_;
    double rmax_;
    double rmin_;
    double z0_;
    double w0_;
    double theta0_;
    double q0_;
    double psi0_;
    double r0_;
    double delta_v0_;
    double delta_s0_;
    double delta_m0_;
    double delta_h0_;
    double zf_;
    double thetaf_;
    double psif_;
    double t0_;
    double tend_;
    double dv_prev_;
    double ds_prev_;
    double dm_prev_;
    double dh_prev_;
    double dv_dot_max_;
    double dm_dot_max_;
    double ds_dot_max_;
    double dh_dot_max_;

    std::array<std::array<double, 4>, 4> A_;
    std::array<std::array<double, 3>, 4> B_;
    std::array<std::array<double, 2>, 2> C_;
    std::array<std::array<double, 1>, 2> D_;
    Bebot bebot_;
    std::vector<Number> solution_u_;
    std::vector<Number> solution_x2_;
    std::vector<Number> solution_x_;
    Number final_obj_value_;
    std::vector<double> final_time_;
    std::vector<std::vector<double>> bernsteinpoly_resultu_;
    std::vector<std::vector<double>> bernsteinpoly_resultx2_;
    std::vector<std::vector<double>> bernsteinpoly_resultz_;

    // ====== CACHES (computed once) ======
    std::vector<double> Dm_cached_;
    std::vector<double> tau_vec_cached_;
    std::vector<double> r_tau_row_cached_;
    std::vector<double> ones_cached_;
    std::vector<double> zf_vector_cached_;
    std::vector<double> thetaf_vector_cached_;
    std::vector<double> psif_vector_cached_;

    void build_caches_() {
        const int Np = N_ + 1;

        // Cache Dm (BeBOT)
        Dm_cached_ = bebot_.getDifferentiationMatrix();

        // Cache tau_vec and r_tau_row (same as your eval_g code)
        tau_vec_cached_.resize(Np);
        r_tau_row_cached_.resize(Np);

        if (Np == 1) {
            tau_vec_cached_[0] = -1.0;
            const double denom = 1.0 - tau_vec_cached_[0];
            r_tau_row_cached_[0] = 2.0 / (denom * denom);
        } else {
            const double tau0 = -1.0;
            const double tau1 = 0.999;
            const double step = (tau1 - tau0) / static_cast<double>(N_);
            for (int i = 0; i < Np; ++i) {
                tau_vec_cached_[i] = tau0 + step * static_cast<double>(i);
                const double denom = 1.0 - tau_vec_cached_[i];
                r_tau_row_cached_[i] = 2.0 / (denom * denom);
            }
        }

        // Cache ones
        ones_cached_.assign(Np, 1.0);

        // Cache target vectors for eval_f
        zf_vector_cached_.assign(Np, zf_);
        thetaf_vector_cached_.assign(Np, thetaf_);
        psif_vector_cached_.assign(Np, psif_);
    }

public:
    const std::vector<std::vector<double>>& get_bernsteinpoly_result() const {
        return bernsteinpoly_resultz_;
    }
};

extern "C" {
    PointSetProblem* create_point_set_problem(int N, double tf, double delta_v_max, double delta_v_min,
        double delta_s_max, double delta_s_min, double delta_m_max, double delta_m_min, double delta_h_max, double delta_h_min,
        double zmax, double zmin, double wmax, double wmin, double thetamax, double thetamin, double qmax, double qmin,
        double psimax, double psimin, double rmax, double rmin,
        double z0, double w0, double theta0, double q0,
        double psi0, double r0,
        double delta_v0, double delta_s0, double delta_m0, double delta_h0,
        double zf, double thetaf, double psif,
        double a11, double a12, double a13, double a14,
        double a21, double a22, double a23, double a24,
        double a31, double a32, double a33, double a34,
        double a41, double a42, double a43, double a44,
        double b11, double b12, double b13,
        double b21, double b22, double b23,
        double b31, double b32, double b33,
        double b41, double b42, double b43,
        double c11, double c12,
        double c21, double c22,
        double d11,
        double d21,
        double t0, double tend,
        double dv_prev, double ds_prev, double dm_prev, double dh_prev,
        double dv_dot_max, double dm_dot_max, double ds_dot_max, double dh_dot_max) {

        return new PointSetProblem(N, tf, delta_v_max, delta_v_min, delta_s_max, delta_s_min,
            delta_m_max, delta_m_min, delta_h_max, delta_h_min,
            zmax, zmin, wmax, wmin, thetamax, thetamin, qmax, qmin,
            psimax, psimin, rmax, rmin,
            z0, w0, theta0, q0, psi0, r0,
            delta_v0, delta_s0, delta_m0, delta_h0,
            zf, thetaf, psif,
            a11, a12, a13, a14,
            a21, a22, a23, a24,
            a31, a32, a33, a34,
            a41, a42, a43, a44,
            b11, b12, b13,
            b21, b22, b23,
            b31, b32, b33,
            b41, b42, b43,
            c11, c12,
            c21, c22,
            d11,
            d21,
            t0, tend, dv_prev, ds_prev, dm_prev, dh_prev, dv_dot_max, dm_dot_max, ds_dot_max, dh_dot_max);
    }

    void solve_point_set_problem(PointSetProblem* problem) {
        SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
        app->Options()->SetStringValue("linear_solver", "ma57");
        app->Options()->SetStringValue("mu_strategy", "adaptive");
        app->Options()->SetStringValue("gradient_approximation", "finite-difference-values");
        app->Options()->SetStringValue("jacobian_approximation", "finite-difference-values");
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");
        app->Options()->SetIntegerValue("max_iter", 400);
        app->Options()->SetNumericValue("tol",             1e-6);
        app->Options()->SetNumericValue("constr_viol_tol", 1e-6);
        app->Options()->SetNumericValue("acceptable_tol",  1e-6);
        app->Options()->SetIntegerValue("print_level", 5);

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
            for (Index i = 0; i < (Index)solution_x.size(); i++) {
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
        std::copy(sol.begin(), sol.end(), solution);
    }

    double get_final_objective_value(PointSetProblem* problem) {
        return problem->get_final_obj_value();
    }

    void destroy_point_set_problem(PointSetProblem* problem) {
        delete problem;
    }
}

// g++ -shared -fPIC -o libbebot_mpc_auv_v1_threed_nonlin_psi.so ~/dev/optimization/BeBOT_cpp_v2/examples/bebot/ma_57/example_mpc_auv_threed/libbebot_mpc_auv_v_nonlin_psi.cpp ~/dev/optimization/BeBOT_cpp_v2/examples/bebot/ma_57/example_mpc_auv_threed/state_space_matrices.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bebot.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinpoly.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteindifferentialmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinmatrix_a2b.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/degelevmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/nchoosek_mod.cpp -I~/dev/optimization/BeBOT_cpp_v2/include -I./Ipopt/src/ -L./Ipopt/src/.libs -lipopt -L/opt/intel/oneapi/mkl/latest/lib/intel64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -ldl -lm -lpthread -lstdc++


// g++ -shared -fPIC -o libbebot_mpc_auv_v1_threed_nonlin_psi.so ~/dev/optimization/BeBOT_cpp_v2/examples/bebot/ma_57/example_mpc_auv_threed/libbebot_mpc_auv_v_nonlin_psi.cpp ~/dev/optimization/BeBOT_cpp_v2/examples/bebot/ma_57/example_mpc_auv_threed/state_space_matrices.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bebot.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinpoly.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteindifferentialmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinmatrix_a2b.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/degelevmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/nchoosek_mod.cpp -I~/dev/optimization/BeBOT_cpp_v2/include -I./Ipopt/src/ -L./Ipopt/src/.libs -lipopt -L/opt/intel/oneapi/mkl/latest/lib/intel64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -ldl -lm -lpthread -lstdc++
