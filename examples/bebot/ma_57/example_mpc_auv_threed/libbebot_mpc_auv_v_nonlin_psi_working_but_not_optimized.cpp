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
        // double delta_n_max, double delta_n_min, 
        double zmax, double zmin, double wmax, double wmin, double thetamax, double thetamin, double qmax, double qmin,
        // double umax, double umin, 
        double psimax, double psimin, double rmax, double rmin,
        // double xmax, double xmin, double ymax, double ymin, 
        double z0, double w0, double theta0, double q0,
        // double u0, 
        double psi0,  double r0,
        // double x0,  double y0, 
        double delta_v0, double delta_s0, double delta_m0, double delta_h0, 
        // double delta_n0,
        double zf, double thetaf, 
        // double xf, double yf, 
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
        // delta_n_max_(delta_n_max), delta_n_min_(delta_n_min), 
        zmax_(zmax), zmin_(zmin), wmax_(wmax), wmin_(wmin), thetamax_(thetamax), thetamin_(thetamin), qmax_(qmax), qmin_(qmin),
        // umax_(umax), umin_(umin), 
        psimax_(psimax), psimin_(psimin), rmax_(rmax), rmin_(rmin), 
        // xmax_(xmax), xmin_(xmin), ymax_(ymax), ymin_(ymin), 
        z0_(z0), w0_(w0), theta0_(theta0), q0_(q0),
        // u0_(u0), 
        psi0_(psi0), r0_(r0),
        // x0_(x0), y0_(y0), 
        delta_v0_(delta_v0), delta_s0_(delta_s0), delta_m0_(delta_m0), delta_h0_(delta_h0), 
        // delta_n0_(delta_n0),
        zf_(zf), thetaf_(thetaf), 
        // xf_(xf), yf_(yf), 
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

        // Construct the B matrix
        D_ = {{
            {d11},
            {d21}
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
        n = 10 * (N_ + 1); // z, theta, w, q, psi, r, dv, dm, ds, dh 
        m = 10 * (N_ + 1); // 6*(N+1) dynamics + 4*(N+1) control derivatives
        nnz_jac_g = n * m;  
        nnz_h_lag = 0; 
        index_style = TNLP::C_STYLE;
        return true;
    }

    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u) {
        // Precompute bounds to avoid redundant calculations
        std::vector<double> x_lower(n, -std::numeric_limits<double>::infinity());
        std::vector<double> x_upper(n, std::numeric_limits<double>::infinity());
        // z
        for (int i = 1; i < N_ + 1; ++i) {
            x_lower[i] = zmin_;
            x_upper[i] = zmax_;
        }
        x_lower[0] = x_upper[0] = z0_;
        //x_lower[N_] = x_upper[N_] = zf_;

        // theta
        for (int i = N_ + 1 + 1; i < 2 * (N_ + 1); ++i) {
            x_lower[i] = thetamin_;
            x_upper[i] = thetamax_;
        }
        x_lower[(N_ + 1)] = x_upper[(N_ + 1)] = theta0_;

        // w
        for (int i = 2 * (N_ + 1) + 1; i < 3 * (N_ + 1); ++i) {
            x_lower[i] = wmin_;
            x_upper[i] = wmax_;
        }
        x_lower[2 * (N_ + 1)] = x_upper[2 * (N_ + 1)] = w0_;

        // for (int i = 2 * (N_ + 1) + 1; i < 3 * (N_ + 1); ++i) {
        //     x_lower[i] = thetamin_;
        //     x_upper[i] = thetamax_;
        // }
        // x_lower[2 * (N_ + 1)] = x_upper[2 * (N_ + 1)] = theta0_;
        
        // q
        for (int i = 3 * (N_ + 1) + 1; i < 4 * (N_ + 1); ++i) {
            x_lower[i] = qmin_;
            x_upper[i] = qmax_;
        }
        x_lower[3 * (N_ + 1)] = x_upper[3 * (N_ + 1)] = q0_;//

        // // 
        // for (int i = 4 * (N_ + 1) + 1; i < 5 * (N_ + 1); ++i) {
        //     x_lower[i] = umin_;
        //     x_upper[i] = umax_;
        // }
        // x_lower[4 * (N_ + 1)] = x_upper[4 * (N_ + 1)] = -u0_;

        // psi
        for (int i = 4 * (N_ + 1) + 1; i < 5 * (N_ + 1); ++i) {
            x_lower[i] = psimin_;
            x_upper[i] = psimax_;
        }
        x_lower[4 * (N_ + 1)] = x_upper[4 * (N_ + 1)] = psi0_;

        // r
        for (int i = 5 * (N_ + 1) + 1; i < 6 * (N_ + 1); ++i) {
            x_lower[i] = rmin_;
            x_upper[i] = rmax_;
        }
        x_lower[5 * (N_ + 1)] = x_upper[5 * (N_ + 1)] = r0_;

        // x
        // for (int i = 7 * (N_ + 1) + 1; i < 8 * (N_ + 1); ++i) {
        //     x_lower[i] = xmin_;
        //     x_upper[i] = xmax_;
        // }
        // x_lower[7 * (N_ + 1)] = x_upper[7 * (N_ + 1)] = x0_;
        //x_lower[7 * (N_ + 1)] = x_upper[7 * (N_ + 1)] = r0_;

        // y
        // for (int i = 8 * (N_ + 1); i < 9 * (N_ + 1); ++i) { // 8 * (N_ + 1) + 1
        //     x_lower[i] = ymin_;
        //     x_upper[i] = ymax_;
        // }
        // x_lower[8 * (N_ + 1)] = x_upper[8 * (N_ + 1)] = y0_;
        // //x_lower[8 * (N_ + 1)] = x_upper[8 * (N_ + 1)] = delta_v0_;


        // control input
        for (int i = 6 * (N_ + 1); i < 7 * (N_ + 1); ++i) { // 9 * (N_ + 1) + 1
            x_lower[i] = delta_v_min_;
            x_upper[i] = delta_v_max_;
        }
        x_lower[6 * (N_ + 1)] = x_upper[6 * (N_ + 1)] = dv_prev_;

        for (int i = 7 * (N_ + 1); i < 8 * (N_ + 1); ++i) { // 10 * (N_ + 1) + 1
            x_lower[i] = delta_m_min_;
            x_upper[i] = delta_m_max_;
        }
        x_lower[7 * (N_ + 1)] = x_upper[7 * (N_ + 1)] = dm_prev_;

        for (int i = 8 * (N_ + 1); i < 9 * (N_ + 1); ++i) { // 11 * (N_ + 1) + 1
            x_lower[i] = delta_s_min_;
            x_upper[i] = delta_s_max_;
        }
        x_lower[8 * (N_ + 1)] = x_upper[8 * (N_ + 1)] = ds_prev_;

        for (int i = 9 * (N_ + 1); i < 10 * (N_ + 1); ++i) { // 11 * (N_ + 1) + 1
            x_lower[i] = delta_h_min_;
            x_upper[i] = delta_h_max_;
        }
        x_lower[9 * (N_ + 1)] = x_upper[9 * (N_ + 1)] = dh_prev_;

        // for (int i = 13 * (N_ + 1); i < 14 * (N_ + 1); ++i) { // 11 * (N_ + 1) + 1
        //     x_lower[i] = delta_n_min_;
        //     x_upper[i] = delta_n_max_;
        // }
        //x_lower[11 * (N_ + 1)] = x_upper[10 * (N_ + 1)] = delta_h0_;

        // x_lower[n - 1] = tf_;
        // x_upper[n - 1] = tf_;
        std::copy(x_lower.begin(), x_lower.end(), x_l);
        std::copy(x_upper.begin(), x_upper.end(), x_u);

        // --------------------
        // g bounds (constraints)
        // m must be 10*(N_+1)
        // --------------------

        // dynamics equalities: 6*(N+1) constraints -> g == 0
        for (int i = 0; i < 6 * (N_ + 1); ++i) {
            g_l[i] = 0.0;
            g_u[i] = 0.0;
        }

        // dv_dt bounds: [-0.60, 0.60]
        for (int i = 6 * (N_ + 1); i < 7 * (N_ + 1); ++i) {
            g_l[i] = -dv_dot_max_;
            g_u[i] =  dv_dot_max_;
        }

        // dm_dt bounds: [-8.5, 8.5]
        for (int i = 7 * (N_ + 1); i < 8 * (N_ + 1); ++i) {
            g_l[i] = -dm_dot_max_;
            g_u[i] =  dm_dot_max_;
        }

        // ds_dt bounds: [-0.60, 0.60]
        for (int i = 8 * (N_ + 1); i < 9 * (N_ + 1); ++i) {
            g_l[i] = -ds_dot_max_;
            g_u[i] =  ds_dot_max_;
        }

        // dh_dt bounds: [-0.60, 0.60]
        for (int i = 9 * (N_ + 1); i < 10 * (N_ + 1); ++i) {
            g_l[i] = -dh_dot_max_;
            g_u[i] =  dh_dot_max_;
        }

        // for (int i = 0; i < N_; ++i) {
        //     g_l[9 * (N_ + 1) + i] = -std::numeric_limits<double>::infinity();
        //     g_u[9 * (N_ + 1) + i] = 0.0;
        // }

        // for (Index i = 0; i < n; ++i) std::cout << "x_l["<<i<<"]="<<x_l[i]<<", x_u["<<i<<"]="<<x_u[i]<<"\n";
        // for (Index i = 0; i < m; ++i) std::cout << "g_l["<<i<<"]="<<g_l[i]<<", g_u["<<i<<"]="<<g_u[i]<<"\n";


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
       
        const double wz = 100;
        const double wtheta = 1;
        const double wpsi = 10;

        const double wdv = 100;
        const double wdm = 0.01;
        const double wds = 0.01;
        const double wdh = 0.001;

        obj_value = 0.0;

        std::vector<double> zf_vector(N_ + 1, zf_);
        std::vector<double> thetaf_vector(N_ + 1, thetaf_);
        std::vector<double> psif_vector(N_ + 1, psif_);
        

        std::vector<double> z_vector(x, x + (N_ + 1));
        std::vector<double> theta_vector(x + 1 * (N_ + 1), x + 2 * (N_ + 1));
        std::vector<double> psi_vector(x + 4 * (N_ + 1), x + 5 * (N_ + 1));

        std::vector<double> dv_vector(x + 6 * (N_ + 1), x + 7 * (N_ + 1));
        std::vector<double> dm_vector(x + 7 * (N_ + 1), x + 8 * (N_ + 1));
        std::vector<double> ds_vector(x + 8 * (N_ + 1), x + 9 * (N_ + 1));
        std::vector<double> dh_vector(x + 9 * (N_ + 1), x + 10 * (N_ + 1));

        // --- errors ---
        std::vector<double> z_diff(N_ + 1);
        std::vector<double> theta_diff(N_ + 1);
        std::vector<double> psi_diff(N_ + 1);

        vdSub(N_ + 1, z_vector.data(),     zf_vector.data(),     z_diff.data());
        vdSub(N_ + 1, theta_vector.data(), thetaf_vector.data(), theta_diff.data());
        vdSub(N_ + 1, psi_vector.data(),   psif_vector.data(),   psi_diff.data());

        // for(int i=0;i<N_+1;++i) std::cout<<"z_vector["<<i<<"]="<<z_vector[i]<<" zf_vector["<<i<<"]="<<zf_vector[i]<<" z_diff["<<i<<"]="<<z_diff[i]<<"\n";
        // for(int i=0;i<N_+1;++i) std::cout<<"theta_vector["<<i<<"]="<<theta_vector[i]<<" thetaf_vector["<<i<<"]="<<thetaf_vector[i]<<" theta_diff["<<i<<"]="<<theta_diff[i]<<"\n";
        // for(int i=0;i<N_+1;++i) std::cout<<"psi_vector["<<i<<"]="<<psi_vector[i]<<" psif_vector["<<i<<"]="<<psif_vector[i]<<" psi_diff["<<i<<"]="<<psi_diff[i]<<"\n";


        // --- squared errors ---
        std::vector<double> z_diff_sqr(N_ + 1);
        std::vector<double> theta_diff_sqr(N_ + 1);
        std::vector<double> psi_diff_sqr(N_ + 1);

        vdSqr(N_ + 1, z_diff.data(),     z_diff_sqr.data());
        vdSqr(N_ + 1, theta_diff.data(), theta_diff_sqr.data());
        vdSqr(N_ + 1, psi_diff.data(),   psi_diff_sqr.data());

        // for(int i=0;i<N_+1;++i) std::cout<<"z_diff["<<i<<"]="<<z_diff[i]<<" z_diff_sqr["<<i<<"]="<<z_diff_sqr[i]<<"\n";
        // for(int i=0;i<N_+1;++i) std::cout<<"theta_diff["<<i<<"]="<<theta_diff[i]<<" theta_diff_sqr["<<i<<"]="<<theta_diff_sqr[i]<<"\n";
        // for(int i=0;i<N_+1;++i) std::cout<<"psi_diff["<<i<<"]="<<psi_diff[i]<<" psi_diff_sqr["<<i<<"]="<<psi_diff_sqr[i]<<"\n";


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
            const double tau1 = 0.999;
            const double step = (tau1 - tau0) / static_cast<double>(N_); // (N+1)-1 = N
            for (int i = 0; i < (N_ + 1); ++i) {
                const double tau = tau0 + step * static_cast<double>(i);
                const double denom = (1.0 - tau);
                r_tau[i] = 2.0 / (denom * denom);
            }
        }
        // for (int i = 0; i < N_ + 1; ++i) std::cout << "r_tau[" << i << "]=" << r_tau[i] << "\n";

        // --- elementwise weighted squared errors: r_tau .* (error.^2) ---
        std::vector<double> z_weighted(N_ + 1);
        std::vector<double> theta_weighted(N_ + 1);
        std::vector<double> psi_weighted(N_ + 1);

        vdMul(N_ + 1, r_tau.data(), z_diff_sqr.data(),     z_weighted.data());
        vdMul(N_ + 1, r_tau.data(), theta_diff_sqr.data(), theta_weighted.data());
        vdMul(N_ + 1, r_tau.data(), psi_diff_sqr.data(),   psi_weighted.data());

        // for (int i = 0; i < N_ + 1; ++i) std::cout << "r_tau[" << i << "]=" << r_tau[i] << (i+1<N_+1?" ":"\n");
        // for (int i = 0; i < N_ + 1; ++i) std::cout << "z_diff_sqr[" << i << "]=" << z_diff_sqr[i] << (i+1<N_+1?" ":"\n");
        // for (int i = 0; i < N_ + 1; ++i) std::cout << "z_weighted[" << i << "]=" << z_weighted[i] << (i+1<N_+1?" ":"\n");

        // for (int i = 0; i < N_ + 1; ++i) std::cout << "theta_diff_sqr[" << i << "]=" << theta_diff_sqr[i] << (i+1<N_+1?" ":"\n");
        // for (int i = 0; i < N_ + 1; ++i) std::cout << "theta_weighted[" << i << "]=" << theta_weighted[i] << (i+1<N_+1?" ":"\n");

        // for (int i = 0; i < N_ + 1; ++i) std::cout << "psi_diff_sqr[" << i << "]=" << psi_diff_sqr[i] << (i+1<N_+1?" ":"\n");
        // for (int i = 0; i < N_ + 1; ++i) std::cout << "psi_weighted[" << i << "]=" << psi_weighted[i] << (i+1<N_+1?" ":"\n");


        // --- sum() via dot with ones ---
        std::vector<double> ones(N_ + 1, 1.0);

        const double sum_z_state     = cblas_ddot(N_ + 1, z_weighted.data(),     1, ones.data(), 1);
        const double sum_theta_state = cblas_ddot(N_ + 1, theta_weighted.data(), 1, ones.data(), 1);
        const double sum_psi_state   = cblas_ddot(N_ + 1, psi_weighted.data(),   1, ones.data(), 1);

        // for (int i = 0; i < N_ + 1; ++i) std::cout << "z_weighted[" << i << "]=" << z_weighted[i] << (i+1<N_+1?" ":"\n");
        // for (int i = 0; i < N_ + 1; ++i) std::cout << "ones[" << i << "]=" << ones[i] << (i+1<N_+1?" ":"\n");
        // std::cout << "sum_z_state=" << sum_z_state << "\n";

        const double state_term =
            wz     * sum_z_state +
            wtheta * sum_theta_state +
            wpsi   * sum_psi_state;

        // --- control term: sum(u.^2) (UNWEIGHTED, same as MATLAB snippet) ---

        std::vector<double> dv_sqr(N_ + 1);
        std::vector<double> dm_sqr(N_ + 1);
        std::vector<double> ds_sqr(N_ + 1);
        std::vector<double> dh_sqr(N_ + 1);

        vdSqr(N_ + 1, dv_vector.data(), dv_sqr.data());
        vdSqr(N_ + 1, dm_vector.data(), dm_sqr.data());
        vdSqr(N_ + 1, ds_vector.data(), ds_sqr.data());
        vdSqr(N_ + 1, dh_vector.data(), dh_sqr.data());


        const double sum_dv = cblas_ddot(N_ + 1, dv_sqr.data(), 1, ones.data(), 1);
        const double sum_dm = cblas_ddot(N_ + 1, dm_sqr.data(), 1, ones.data(), 1);
        const double sum_ds = cblas_ddot(N_ + 1, ds_sqr.data(), 1, ones.data(), 1);
        const double sum_dh = cblas_ddot(N_ + 1, dh_sqr.data(), 1, ones.data(), 1);

        const double control_term =
            wdv * sum_dv +
            wdm * sum_dm +
            wds * sum_ds +
            wdh * sum_dh;

        // --- DEBUG dv: dv_vector -> dv_sqr -> ddot with ones -> contribution ---
        // for (int i = 0; i < N_ + 1; ++i) std::cout << "dv_vector[" << i << "]=" << dv_vector[i] << (i+1<N_+1?" ":"\n");
        // for (int i = 0; i < N_ + 1; ++i) std::cout << "dv_sqr["    << i << "]=" << dv_sqr[i]    << (i+1<N_+1?" ":"\n");
        // for (int i = 0; i < N_ + 1; ++i) std::cout << "ones["      << i << "]=" << ones[i]      << (i+1<N_+1?" ":"\n");
        // std::cout << "sum_dv=" << sum_dv << "\n";
        // std::cout << "wdv=" << wdv << "\n";
        // std::cout << "wdv*sum_dv=" << (wdv * sum_dv) << "\n";


        obj_value = state_term + control_term;

        return true;
    }

    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
        
        // eval_g implements constraints g(x):
        //   - 6*(N+1) dynamics equalities (BeBOT tau-derivative form)
        //   - 4*(N+1) control time-derivative bounds (dv_dt, dm_dt, ds_dt, dh_dt)
        //
        // MATLAB reference:
        //   Xd_tau = Xd * Dm
        //   Xh_tau = Xh * Dm
        //   res_depth = Xd_tau - (A*Xd + B*Ud) .* r_tauRow
        //   res_horiz = Xh_tau - (C*Xh + D*Uh) .* r_tauRow
        //   u_tau = u_row * Dm
        //   u_dt  = u_tau ./ r_tauRow
        //
        // Your IPOPT bounds setup:
        //   g[0..6*(N+1)-1] == 0
        //   g[6*(N+1)..10*(N+1)-1] in [-u_dot_max, +u_dot_max]
        // ============================================================

        // -----------------------------
        // 0) Get BeBOT differentiation matrix Dm for (N_, tf_)
        // MATLAB: [~,~,Dm_cache] = BeBOT(N, tf)
        
        Bebot Bebot(N_, tf_);
        Bebot.calculate();
        const auto& Dm = Bebot.getDifferentiationMatrix();
        
        //std::cout << "tf (x[n-1]) = " << x[n - 1] << std::endl;

        // Assume Dm is a flat vector of length (N_+1)*(N_+1).  Print it as an (N_+1)x(N_+1) matrix:
        // std::cout << "Dm = " << std::endl;
        // for (int i = 0; i < N_ + 1; ++i) {
        //     for (int j = 0; j < N_ + 1; ++j) {
        //         std::cout << Dm[i * (N_ + 1) + j] << "  ";
        //     }
        //     std::cout << std::endl;
        // }


        // -----------------------------
        // 1) Build tau and r_tau_row (same as MATLAB else-branch)
        // MATLAB:
        //   tau   = linspace(-1, 0.999, N+1)
        //   r_tau = (2 ./ (1 - tau).^2).'
        // -----------------------------
        std::vector<double> tau_vec(N_ + 1);
        std::vector<double> r_tau_row(N_ + 1);

        if (N_ + 1 == 1) {
            tau_vec[0] = -1.0;
            const double denom = 1.0 - tau_vec[0];
            r_tau_row[0] = 2.0 / (denom * denom);
        } else {
            const double tau0 = -1.0;
            const double tau1 = 0.999;
            const double step = (tau1 - tau0) / static_cast<double>(N_); // (N+1)-1 = N
            for (int i = 0; i < (N_ + 1); ++i) {
                tau_vec[i] = tau0 + step * static_cast<double>(i);
                const double denom = 1.0 - tau_vec[i];
                r_tau_row[i] = 2.0 / (denom * denom);
            }
        }

        // -----------------------------
        // 2) Unpack decision vector x into state/control knot vectors
        // X = [z, theta, w, q, psi, r, dv, dm, ds, dh] each length (N+1)
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
        // Your convention here (using Dm^T in gemv) should match your Dm storage.
        // dyn1..dyn6 represent z_tau, theta_tau, w_tau, q_tau, psi_tau, r_tau.
        // -----------------------------
        std::vector<double> dyn1(N_ + 1);
        std::vector<double> dyn2(N_ + 1);
        std::vector<double> dyn3(N_ + 1);
        std::vector<double> dyn4(N_ + 1);
        std::vector<double> dyn5(N_ + 1);
        std::vector<double> dyn6(N_ + 1);
        //std::cout << "z_vector = ["; for (int i=0;i<N_+1;++i) std::cout << z_vector[i] << (i<N_? ", ":""); std::cout << "]\n";
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, z_vector.data(), 1, 0.0, dyn1.data(), 1); // z_tau
        //std::cout << "dyn1 = ["; for (int i=0;i<N_+1;++i) std::cout << dyn1[i] << (i<N_? ", ":""); std::cout << "]\n";
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, theta_vector.data(), 1, 0.0, dyn2.data(), 1); // theta_tau
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, w_vector.data(), 1, 0.0, dyn3.data(), 1); // w_tau
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, q_vector.data(), 1, 0.0, dyn4.data(), 1); // q_tau
        
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, psi_vector.data(), 1, 0.0, dyn5.data(), 1); // psi_tau
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, r_vector.data(), 1, 0.0, dyn6.data(), 1); // r_tau
        

        // -----------------------------
        // 4) Build depth subsystem matrices Xd (4x(N+1)) and Ud (3x(N+1))
        // Xd = [z; theta; w; q], Ud = [dv; dm; ds]
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
        // MATLAB: Ad_cache*Xd + Bd_cache*Ud  -> 4x(N+1)
        // We store it as X2_matrix_flat (4 rows, each length (N+1)).
        // -----------------------------
        std::vector<double> X2_matrix_flat(4 * (N_ + 1), 0.0);

        // X2 = A*Xd
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4, (N_ + 1), 4, 1.0, &A_[0][0], 4, X1_matrix_flat.data(), (N_ + 1), 0.0, X2_matrix_flat.data(), (N_ + 1));

        // X2 += B*Ud
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4, (N_ + 1), 3, 1.0, &B_[0][0], 3,U_matrix_flat.data(), (N_ + 1),1.0, X2_matrix_flat.data(), (N_ + 1));


        // auto pv = [&](const char* name, const std::vector<double>& v){
        //     std::cout << name << " = [";
        //     for (int i=0;i<(int)v.size();++i) std::cout << v[i] << (i+1<(int)v.size()? ", ":"");
        //     std::cout << "]\n";
        // };

        // auto pm = [&](const char* name, const double* M, int R, int C, bool rowMajor=true){
        //     std::cout << name << " ("<<R<<"x"<<C<<")\n";
        //     for(int r=0;r<R;++r){
        //         for(int c=0;c<C;++c){
        //             const int idx = rowMajor ? (r*C + c) : (c*R + r);
        //             std::cout << M[idx] << (c+1<C? "  ":"");
        //         }
        //         std::cout << "\n";
        //     }
        // };

        // pm("A_", &A_[0][0], 4, 4, /*rowMajor=*/true);
        // pm("B_", &B_[0][0], 4, 3, /*rowMajor=*/true);

        // pm("Xd (X1_matrix_flat)", X1_matrix_flat.data(), 4, (N_+1), /*rowMajor=*/true);
        // pm("Ud (U_matrix_flat)",  U_matrix_flat.data(),  3, (N_+1), /*rowMajor=*/true);

        // /* --- after first dgemm (X2 = A*Xd) --- */
        // pm("X2 after A*Xd", X2_matrix_flat.data(), 4, (N_+1), /*rowMajor=*/true);

        // /* --- after second dgemm (X2 += B*Ud) --- */
        // pm("X2 after +B*Ud", X2_matrix_flat.data(), 4, (N_+1), /*rowMajor=*/true);


        // -----------------------------
        // 6) Build horizontal subsystem matrices Xh (2x(N+1)) and Uh (1x(N+1))
        // Xh = [psi; r], Uh = [dh]
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
        // MATLAB: Ah_cache*Xh + Bh_cache*Uh -> 2x(N+1)
        // Here C_ is 2x2, D_ is 2x1, so result is 2x(N+1).
        // -----------------------------
        std::vector<double> X2_matrix_flat_cd(2 * (N_ + 1), 0.0);

        // // ---- prints for horizontal subsystem ----
        // pm("C_", &C_[0][0], 2, 2, /*rowMajor=*/true);
        // pm("D_", &D_[0][0], 2, 1, /*rowMajor=*/true);

        // pm("Xh (X1_matrix_flat_cd)", X1_matrix_flat_cd.data(), 2, (N_+1), /*rowMajor=*/true);
        // pm("Uh (U_matrix_flat_cd)",  U_matrix_flat_cd.data(),  1, (N_+1), /*rowMajor=*/true);

        // X2_cd = C*Xh
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, (N_ + 1), 2, 1.0, &C_[0][0], 2, X1_matrix_flat_cd.data(), (N_ + 1),0.0, X2_matrix_flat_cd.data(), (N_ + 1));
        // pm("X2_cd after C*Xh", X2_matrix_flat_cd.data(), 2, (N_+1), /*rowMajor=*/true);
        // X2_cd += D*Uh
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, (N_ + 1), 1, 1.0, &D_[0][0], 1, U_matrix_flat_cd.data(), (N_ + 1),1.0, X2_matrix_flat_cd.data(), (N_ + 1));
        // pm("X2_cd after +D*Uh", X2_matrix_flat_cd.data(), 2, (N_+1), /*rowMajor=*/true);
        

        // -----------------------------
        // 8) Unpack RHS rows from X2 matrices into per-state vectors
        // Depth RHS rows: z_rhs, theta_rhs, w_rhs, q_rhs
        // Horizontal RHS rows: psi_rhs, r_rhs
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
        // 9) Apply tau-weighting to RHS: (A*X + B*U) .* r_tau_row
        // MATLAB broadcast:
        //   each column k is multiplied by r_tauRow(k)
        // Here it is elementwise vector multiply per-row.
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

        // // ---- print raw RHS matrices ----
        // pm("X2_matrix_flat (depth RHS) [4x(N+1)]",    X2_matrix_flat.data(),    4, (N_+1), /*rowMajor=*/true);
        // pm("X2_matrix_flat_cd (horiz RHS) [2x(N+1)]", X2_matrix_flat_cd.data(), 2, (N_+1), /*rowMajor=*/true);

        // // ---- extracted RHS vectors ----
        // pv("z_rhs",     z_rhs);
        // pv("theta_rhs", theta_rhs);
        // pv("w_rhs",     w_rhs);
        // pv("q_rhs",     q_rhs);
        // pv("psi_rhs",   psi_rhs);
        // pv("r_rhs",     r_rhs);

        // ---- tau weighting vector ----
        // pv("r_tau_row", r_tau_row);

        // ---- scaled RHS vectors ----
        // pv("z_rhs_scaled",     z_rhs_scaled);
        // pv("theta_rhs_scaled", theta_rhs_scaled);
        // pv("w_rhs_scaled",     w_rhs_scaled);
        // pv("q_rhs_scaled",     q_rhs_scaled);
        // pv("psi_rhs_scaled",   psi_rhs_scaled);
        // pv("r_rhs_scaled",     r_rhs_scaled);


        // -----------------------------
        // 10) Dynamics equality residuals:
        //   res = X_tau - RHS_scaled
        // These must be == 0 in IPOPT bounds.
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

        // ---- inputs to residual vdSub: dyn* and *_rhs_scaled ----
        // pv("dyn1 (z_tau)",     dyn1);
        // pv("z_rhs_scaled",     z_rhs_scaled);
        // pv("res_z",            res_z);

        // pv("dyn2 (theta_tau)", dyn2);
        // pv("theta_rhs_scaled", theta_rhs_scaled);
        // pv("res_theta",        res_theta);

        // pv("dyn3 (w_tau)",     dyn3);
        // pv("w_rhs_scaled",     w_rhs_scaled);
        // pv("res_w",            res_w);

        // pv("dyn4 (q_tau)",     dyn4);
        // pv("q_rhs_scaled",     q_rhs_scaled);
        // pv("res_q",            res_q);

        // pv("dyn5 (psi_tau)",   dyn5);
        // pv("psi_rhs_scaled",   psi_rhs_scaled);
        // pv("res_psi",          res_psi);

        // pv("dyn6 (r_tau)",     dyn6);
        // pv("r_rhs_scaled",     r_rhs_scaled);
        // pv("res_r",            res_r);


        // -----------------------------
        // 11) Control derivatives:
        // MATLAB:
        //   u_tau = Urow * Dm
        //   u_dt  = u_tau ./ r_tauRow
        // We put u_dt directly into dyn7..dyn10 to match your g bounds blocks.
        // -----------------------------
        std::vector<double> dyn7(N_ + 1);   // dv_dt
        std::vector<double> dyn8(N_ + 1);   // dm_dt
        std::vector<double> dyn9(N_ + 1);   // ds_dt
        std::vector<double> dyn10(N_ + 1);  // dh_dt

        std::vector<double> dv_tau(N_ + 1);
        std::vector<double> dm_tau(N_ + 1);
        std::vector<double> ds_tau(N_ + 1);
        std::vector<double> dh_tau(N_ + 1);

        // u_tau = Dm^T * u (consistent with your earlier gemv usage)
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, delta_v_vector.data(), 1, 0.0, dv_tau.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, delta_m_vector.data(), 1, 0.0, dm_tau.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, delta_s_vector.data(), 1, 0.0, ds_tau.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, delta_h_vector.data(), 1, 0.0, dh_tau.data(), 1);

        // u_dt = u_tau ./ r_tau_row
        vdDiv(N_ + 1, dv_tau.data(), r_tau_row.data(), dyn7.data());
        vdDiv(N_ + 1, dm_tau.data(), r_tau_row.data(), dyn8.data());
        vdDiv(N_ + 1, ds_tau.data(), r_tau_row.data(), dyn9.data());
        vdDiv(N_ + 1, dh_tau.data(), r_tau_row.data(), dyn10.data());


        // ---- u_tau = Dm^T * u : print inputs (u) and outputs (u_tau) ----
        // pv("delta_v_vector", delta_v_vector);
        // pv("dv_tau",         dv_tau);

        // pv("delta_m_vector", delta_m_vector);
        // pv("dm_tau",         dm_tau);

        // pv("delta_s_vector", delta_s_vector);
        // pv("ds_tau",         ds_tau);

        // pv("delta_h_vector", delta_h_vector);
        // pv("dh_tau",         dh_tau);

        // // ---- u_dt = u_tau ./ r_tau_row : print r_tau_row and dyn7..dyn10 ----
        // pv("r_tau_row", r_tau_row);

        // pv("dyn7 (dv_dt)",  dyn7);
        // pv("dyn8 (dm_dt)",  dyn8);
        // pv("dyn9 (ds_dt)",  dyn9);
        // pv("dyn10 (dh_dt)", dyn10);

        // -----------------------------
        // 12) Pack constraints into g in the exact block layout:
        //   block 0..5 : dynamics equalities (== 0)
        //   block 6..9 : control derivatives in time (bounded by g_l/g_u)
        // -----------------------------
        for (Index i = 0; i < (N_ + 1); ++i) {
            // dynamics equalities
            g[0 * (N_ + 1) + i] = res_z[i];
            g[1 * (N_ + 1) + i] = res_theta[i];
            g[2 * (N_ + 1) + i] = res_w[i];
            g[3 * (N_ + 1) + i] = res_q[i];
            g[4 * (N_ + 1) + i] = res_psi[i];
            g[5 * (N_ + 1) + i] = res_r[i];

            // control derivative constraints (time derivatives)
            g[6 * (N_ + 1) + i] = dyn7[i];   // dv_dt
            g[7 * (N_ + 1) + i] = dyn8[i];   // dm_dt
            g[8 * (N_ + 1) + i] = dyn9[i];   // ds_dt
            g[9 * (N_ + 1) + i] = dyn10[i];  // dh_dt
        }

        // for (int i = 0; i < N_; ++i) {
        //     g[9 * (N_ + 1) + i] = x_vector[i + 1] - x_vector[i];
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
        solution_x_.resize(10 * (N_ + 1));
        
        // Copy the 10 * (N + 1) elements from the x array
        for (Index i = 0; i < 10 * (N_ + 1); ++i) {
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

        // Helper to print any std::vector<double>
        auto print_vec = [&](const std::string& name, const std::vector<double>& v){
            std::cout << name << " = [";
            for (size_t i = 0; i < v.size(); ++i) {
                std::cout << v[i];
                if (i + 1 < v.size()) std::cout << ", ";
            }
            std::cout << "]\n";
        };

        // After unpacking:
        // print_vec("z_vector",        z_vector);
        // print_vec("theta_vector",    theta_vector);
        // print_vec("w_vector",        w_vector);
        // print_vec("q_vector",        q_vector);
        // print_vec("psi_vector",      psi_vector);
        // print_vec("r_vector",        r_vector);

        // print_vec("delta_v_vector",  delta_v_vector);
        // print_vec("delta_m_vector",  delta_m_vector);
        // print_vec("delta_s_vector",  delta_s_vector);
        // print_vec("delta_h_vector",  delta_h_vector);


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

        // std::vector<std::vector<double>> bernstein_z = BernsteinPoly(z_2d, final_time_, 0, tf_);
        // std::vector<std::vector<double>> bernstein_theta = BernsteinPoly(theta_2d, final_time_, 0, tf_);
        // std::vector<std::vector<double>> bernstein_w = BernsteinPoly(w_2d, final_time_, 0, tf_);
        // std::vector<std::vector<double>> bernstein_q = BernsteinPoly(q_2d, final_time_, 0, tf_);

        // std::vector<std::vector<double>> bernstein_psi = BernsteinPoly(psi_2d, final_time_, 0, tf_);
        // std::vector<std::vector<double>> bernstein_r = BernsteinPoly(r_2d, final_time_, 0, tf_);

        // std::vector<std::vector<double>> bernstein_delta_v = BernsteinPoly(delta_v_2d, final_time_, 0, tf_);
        // std::vector<std::vector<double>> bernstein_delta_m = BernsteinPoly(delta_m_2d, final_time_, 0, tf_);
        // std::vector<std::vector<double>> bernstein_delta_s = BernsteinPoly(delta_s_2d, final_time_, 0, tf_);
        // std::vector<std::vector<double>> bernstein_delta_h = BernsteinPoly(delta_h_2d, final_time_, 0, tf_);

        // 2) Evaluate Bernstein on normalized domain [0,1]
        auto z_norm     = BernsteinPoly(z_2d,     t_norm, 0.0, 1.0);
        auto theta_norm = BernsteinPoly(theta_2d, t_norm, 0.0, 1.0);
        auto w_norm     = BernsteinPoly(w_2d,     t_norm, 0.0, 1.0);
        auto q_norm     = BernsteinPoly(q_2d,     t_norm, 0.0, 1.0);

        auto psi_norm   = BernsteinPoly(psi_2d,   t_norm, 0.0, 1.0);
        auto r_norm     = BernsteinPoly(r_2d,     t_norm, 0.0, 1.0);

        auto dv_norm    = BernsteinPoly(delta_v_2d,    t_norm, 0.0, 1.0);
        auto dm_norm    = BernsteinPoly(delta_m_2d,    t_norm, 0.0, 1.0);
        auto ds_norm    = BernsteinPoly(delta_s_2d,    t_norm, 0.0, 1.0);
        auto dh_norm    = BernsteinPoly(delta_h_2d,    t_norm, 0.0, 1.0);

        // 3) Evaluate Bernstein on real-time domain [0,tf_]
        auto z_real     = BernsteinPoly(z_2d,     t_real, 0.0, tf_);
        auto theta_real = BernsteinPoly(theta_2d, t_real, 0.0, tf_);
        auto w_real     = BernsteinPoly(w_2d,     t_real, 0.0, tf_);
        auto q_real     = BernsteinPoly(q_2d,     t_real, 0.0, tf_);

        auto psi_real   = BernsteinPoly(psi_2d,   t_real, 0.0, tf_);
        auto r_real     = BernsteinPoly(r_2d,     t_real, 0.0, tf_);

        auto dv_real    = BernsteinPoly(delta_v_2d,    t_real, 0.0, tf_);
        auto dm_real    = BernsteinPoly(delta_m_2d,    t_real, 0.0, tf_);
        auto ds_real    = BernsteinPoly(delta_s_2d,    t_real, 0.0, tf_);
        auto dh_real    = BernsteinPoly(delta_h_2d,    t_real, 0.0, tf_);

        auto flatten = [](const std::vector<std::vector<double>>& input) {
            std::vector<double> output;
            for (const auto& row : input) {
                output.insert(output.end(), row.begin(), row.end());
            }
            return output;
        };
        // writeToCSV(final_time_, flatten(bernstein_z), "z.csv");
        // writeToCSV(bebot_.getNodes(), z_vector, "z_controlpoints.csv");
        // writeToCSV(final_time_, flatten(bernstein_w), "w.csv");
        // writeToCSV(bebot_.getNodes(), w_vector, "w_controlpoints.csv");
        // writeToCSV(final_time_, flatten(bernstein_theta), "theta.csv");
        // writeToCSV(bebot_.getNodes(), theta_vector, "theta_controlpoints.csv");
        // writeToCSV(final_time_, flatten(bernstein_q), "q.csv");
        // writeToCSV(bebot_.getNodes(), q_vector, "q_controlpoints.csv");

        // writeToCSV(final_time_, flatten(bernstein_psi), "psi.csv");
        // writeToCSV(bebot_.getNodes(), psi_vector, "psi_controlpoints.csv");
        // writeToCSV(final_time_, flatten(bernstein_r), "r.csv");
        // writeToCSV(bebot_.getNodes(), r_vector, "r_controlpoints.csv");

        // writeToCSV(final_time_, flatten(bernstein_delta_v), "delta_v.csv");
        // writeToCSV(bebot_.getNodes(), delta_v_vector, "delta_v_controlpoints.csv");
        // writeToCSV(final_time_, flatten(bernstein_delta_m), "delta_m.csv");
        // writeToCSV(bebot_.getNodes(), delta_m_vector, "delta_m_controlpoints.csv");
        // writeToCSV(final_time_, flatten(bernstein_delta_s), "delta_s.csv");
        // writeToCSV(bebot_.getNodes(), delta_s_vector, "delta_s_controlpoints.csv");
        // writeToCSV(final_time_, flatten(bernstein_delta_h), "delta_h.csv");
        // writeToCSV(bebot_.getNodes(), delta_h_vector, "delta_h_controlpoints.csv");

        // ---- normalized time outputs (0..1) ----
        // writeToCSV(t_norm, flatten(z_norm),     "z_norm.csv");
        // writeToCSV(bebot_.getNodes(), z_vector, "z_controlpoints_tau.csv");

        // writeToCSV(t_norm, flatten(w_norm),     "w_norm.csv");
        // writeToCSV(bebot_.getNodes(), w_vector, "w_controlpoints_tau.csv");

        // writeToCSV(t_norm, flatten(theta_norm),     "theta_norm.csv");
        // writeToCSV(bebot_.getNodes(), theta_vector, "theta_controlpoints_tau.csv");

        // writeToCSV(t_norm, flatten(q_norm),     "q_norm.csv");
        // writeToCSV(bebot_.getNodes(), q_vector, "q_controlpoints_tau.csv");

        // writeToCSV(t_norm, flatten(psi_norm),     "psi_norm.csv");
        // writeToCSV(bebot_.getNodes(), psi_vector, "psi_controlpoints_tau.csv");

        // writeToCSV(t_norm, flatten(r_norm),     "r_norm.csv");
        // writeToCSV(bebot_.getNodes(), r_vector, "r_controlpoints_tau.csv");

        // writeToCSV(t_norm, flatten(dv_norm),          "delta_v_norm.csv");
        // writeToCSV(bebot_.getNodes(), delta_v_vector, "delta_v_controlpoints_tau.csv");

        // writeToCSV(t_norm, flatten(dm_norm),          "delta_m_norm.csv");
        // writeToCSV(bebot_.getNodes(), delta_m_vector, "delta_m_controlpoints_tau.csv");

        // writeToCSV(t_norm, flatten(ds_norm),          "delta_s_norm.csv");
        // writeToCSV(bebot_.getNodes(), delta_s_vector, "delta_s_controlpoints_tau.csv");

        // writeToCSV(t_norm, flatten(dh_norm),          "delta_h_norm.csv");
        // writeToCSV(bebot_.getNodes(), delta_h_vector, "delta_h_controlpoints_tau.csv");


        // ---- real time outputs (0..tf_) ----
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

        writeToCSV(t_real, flatten(dv_real),         "delta_v_real.csv");
        writeToCSV(bebot_.getNodes(), delta_v_vector,"delta_v_controlpoints_real.csv");

        writeToCSV(t_real, flatten(dm_real),         "delta_m_real.csv");
        writeToCSV(bebot_.getNodes(), delta_m_vector,"delta_m_controlpoints_real.csv");

        writeToCSV(t_real, flatten(ds_real),         "delta_s_real.csv");
        writeToCSV(bebot_.getNodes(), delta_s_vector,"delta_s_controlpoints_real.csv");

        writeToCSV(t_real, flatten(dh_real),         "delta_h_real.csv");
        writeToCSV(bebot_.getNodes(), delta_h_vector,"delta_h_controlpoints_real.csv");


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

public:
    const std::vector<std::vector<double>>& get_bernsteinpoly_result() const { 
        return bernsteinpoly_resultz_; }
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
        
        // std::cout << std::fixed << std::setprecision(5);
        // //Print out every incoming argument:
        // std::cout << "==== create_point_set_problem called ====\n";
        // std::cout << "N = " << N << "\n";
        // std::cout << "tf = " << tf << "\n";

        // std::cout << "delta_v_max = " << delta_v_max
        //           << ", delta_v_min = " << delta_v_min << "\n";
        // std::cout << "delta_s_max = " << delta_s_max
        //           << ", delta_s_min = " << delta_s_min << "\n";
        // std::cout << "delta_m_max = " << delta_m_max
        //           << ", delta_m_min = " << delta_m_min << "\n";
        // std::cout << "delta_h_max = " << delta_h_max
        //           << ", delta_h_min = " << delta_h_min << "\n";


        // std::cout << "zmax = " << zmax << ", zmin = " << zmin << "\n";
        // std::cout << "wmax = " << wmax << ", wmin = " << wmin << "\n";
        // std::cout << "thetamax = " << thetamax << ", thetamin = " << thetamin << "\n";
        // std::cout << "qmax = " << qmax << ", qmin = " << qmin << "\n\n";

        // std::cout << "psimax = " << psimax << ", psimin = " << psimin << "\n";
        // std::cout << "rmax = " << rmax << ", rmin = " << rmin << "\n\n";

        // std::cout << "Initial states:\n";
        // std::cout << "  z0 = " << z0 << ", w0 = " << w0
        //           << ", theta0 = " << theta0 << ", q0 = " << q0 << "\n";
        // std::cout << " psi0 = " << psi0 << ", r0 = " << r0 << "\n";

        // std::cout << "Initial controls:\n";
        // std::cout << "  delta_v0 = " << delta_v0
        //           << ", delta_s0 = " << delta_s0
        //           << ", delta_m0 = " << delta_m0
        //           << ", delta_h0 = " << delta_h0 << "\n\n";

        // std::cout << "Final targets:\n";
        // std::cout << "  zf = " << zf
        //           << ", thetaf = " << thetaf
        //           << ", psif = " << psif << "\n\n";

        // std::cout << "A‐matrix (4×4):\n";
        // std::cout << "  [" << a11 << ", " << a12 << ", " << a13 << ", " << a14 << "]\n";
        // std::cout << "  [" << a21 << ", " << a22 << ", " << a23 << ", " << a24 << "]\n";
        // std::cout << "  [" << a31 << ", " << a32 << ", " << a33 << ", " << a34 << "]\n";
        // std::cout << "  [" << a41 << ", " << a42 << ", " << a43 << ", " << a44 << "]\n\n";

        // std::cout << "B‐matrix (4×3):\n";
        // std::cout << "  [" << b11 << ", " << b12 << ", " << b13 << "]\n";
        // std::cout << "  [" << b21 << ", " << b22 << ", " << b23 << "]\n";
        // std::cout << "  [" << b31 << ", " << b32 << ", " << b33 << "]\n";
        // std::cout << "  [" << b41 << ", " << b42 << ", " << b43 << "]\n\n";

        // std::cout << "C‐matrix (2×2):\n";
        // std::cout << "  [" << c11 << ", " << c12 << "]\n";
        // std::cout << "  [" << c21 << ", " << c22 << "]\n\n";

        // std::cout << "D‐matrix (2×1):\n";
        // std::cout << "  [" << d11 << "]\n";
        // std::cout << "  [" << d21 << "]\n\n";

        // std::cout << "t0 = " << t0 << ", tend = " << tend << "\n";
        // std::cout << "===========================================\n\n";
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
        app->Options()->SetIntegerValue("max_iter", 1000);
        app->Options()->SetNumericValue("tol",             1e-6);   // OptimalityTolerance = 1e-3    1
        app->Options()->SetNumericValue("constr_viol_tol", 1e-6);
        app->Options()->SetNumericValue("acceptable_tol",        1e-6); //1
        // somewhere before app->Initialize():
        //app->Options()->SetIntegerValue("max_line_search_step_retries", 200);
        //app->Options()->SetNumericValue("alpha_for_y", 0.6);
        //app->Options()->SetNumericValue("beta_for_y",  0.4);


        //app->Options()->SetNumericValue("constr_viol_tol", 1e-6);
        app->Options()->SetIntegerValue("print_level", 3); 
        //app->Options()->SetStringValue("nlp_scaling_method", "gradient-based");
        //app->Options()->SetIntegerValue("max_line_search_step_retries", 50);
        //app->Options()->SetNumericValue("alpha_for_y", 0.6);



        //app->Options()->SetNumericValue("finite_difference_rel_step", 1e-4);
        //app->Options()->SetNumericValue("finite_difference_abs_step", 1e-8);
        //app->Options()->SetNumericValue("tol",              1e-3);
        //app->Options()->SetNumericValue("constr_viol_tol",  1e-4);
        //app->Options()->SetNumericValue("dual_inf_tol",     1e-4);

        //---- give the line-search more chances ----
        //app->Options()->SetIntegerValue("max_line_search_step_retries",  20);
        //app->Options()->SetNumericValue("alpha_for_y", 0.6);    // trial fraction for filter
        //app->Options()->SetNumericValue("beta_for_y",  0.4);

        //---- (optionally) loosen “acceptable” termination ----
        //app->Options()->SetNumericValue("acceptable_obj_change_tol", 1e-2);
        //app->Options()->SetIntegerValue("acceptable_iter",      5);

        //app->Options()->SetStringValue ("derivative_test",      "first-order");
        //app->Options()->SetStringValue ("derivative_test_print_all","yes");

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

// g++ -shared -fPIC -o libbebot_mpc_auv_v1_threed_nonlin_psi.so ~/dev/optimization/BeBOT_cpp_v2/examples/bebot/ma_57/example_mpc_auv_threed/libbebot_mpc_auv_v_nonlin_psi.cpp ~/dev/optimization/BeBOT_cpp_v2/examples/bebot/ma_57/example_mpc_auv_threed/state_space_matrices.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bebot.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinpoly.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteindifferentialmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinmatrix_a2b.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/degelevmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/nchoosek_mod.cpp -I~/dev/optimization/BeBOT_cpp_v2/include -I./Ipopt/src/ -L./Ipopt/src/.libs -lipopt -L/opt/intel/oneapi/mkl/latest/lib/intel64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -ldl -lm -lpthread -lstdc++
