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

using namespace Ipopt;

class PointSetProblem : public Ipopt::TNLP {
public:
    PointSetProblem(int N, double tf, double delta_v_max, double delta_v_min, 
        double delta_s_max, double delta_s_min, double delta_m_max, double delta_m_min, double delta_h_max, double delta_h_min, double delta_n_max, double delta_n_min, 
        double zmax, double zmin, double wmax, double wmin, double thetamax, double thetamin, double qmax, double qmin,
        double umax, double umin, double psimax, double psimin, double rmax, double rmin,
        double xmax, double xmin, double ymax, double ymin, 
        double z0, double w0, double theta0, double q0,
        double u0, double psi0,  double r0,
        double x0,  double y0, 
        double delta_v0, double delta_s0, double delta_m0, double delta_h0, double delta_n0,
        double zf, double thetaf, double xf, double yf, double psif,
        double a11, double a12, double a13, double a14, 
        double a21, double a22, double a23, double a24,
        double a31, double a32, double a33, double a34, 
        double a41, double a42, double a43, double a44,  
        double b11, double b12, double b13, 
        double b21, double b22, double b23, 
        double b31, double b32, double b33, 
        double b41, double b42, double b43, 
        double c11, double c12, double c13,  
        double c21, double c22, double c23, 
        double c31, double c32, double c33,  
        double d11, double d12, 
        double d21, double d22, 
        double d31, double d32, 
        double t0, double tend)
        : N_(N), tf_(tf), delta_v_max_(delta_v_max), delta_v_min_(delta_v_min), delta_m_max_(delta_m_max), delta_m_min_(delta_m_min), 
        delta_s_max_(delta_s_max), delta_s_min_(delta_s_min), delta_h_max_(delta_h_max), delta_h_min_(delta_h_min), delta_n_max_(delta_n_max), delta_n_min_(delta_n_min), 
        zmax_(zmax), zmin_(zmin), wmax_(wmax), wmin_(wmin), thetamax_(thetamax), thetamin_(thetamin), qmax_(qmax), qmin_(qmin),
        umax_(umax), umin_(umin), psimax_(psimax), psimin_(psimin), rmax_(rmax), rmin_(rmin), 
        xmax_(xmax), xmin_(xmin), ymax_(ymax), ymin_(ymin), 
        z0_(z0), w0_(w0), theta0_(theta0), q0_(q0),
        u0_(u0), psi0_(psi0), r0_(r0),
        x0_(x0), y0_(y0), 
        delta_v0_(delta_v0), delta_s0_(delta_s0), delta_m0_(delta_m0), delta_h0_(delta_h0), delta_n0_(delta_n0),
        zf_(zf), thetaf_(thetaf), xf_(xf), yf_(yf), psif_(psif), bebot_(N, tf_), t0_(t0), tend_(tend) {
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
            {c11, c12, c13},
            {c21, c22, c23},
            {c31, c32, c33}
        }};

        // Construct the B matrix
        D_ = {{
            {d11, d12},
            {d21, d22},
            {d31, d32}
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
        n = 14 * (N_ + 1) + 1; 
        m = 9 * (N_ + 1);// + N_;
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
        for (int i = 1; i < N_; ++i) {
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

        // y
        for (int i = 4 * (N_ + 1) + 1; i < 5 * (N_ + 1); ++i) {
            x_lower[i] = umin_;
            x_upper[i] = umax_;
        }
        x_lower[4 * (N_ + 1)] = x_upper[4 * (N_ + 1)] = -u0_;

        // psi
        for (int i = 5 * (N_ + 1) + 1; i < 6 * (N_ + 1); ++i) {
            x_lower[i] = psimin_;
            x_upper[i] = psimax_;
        }
        x_lower[5 * (N_ + 1)] = x_upper[5 * (N_ + 1)] = psi0_;

        // v
        for (int i = 6 * (N_ + 1) + 1; i < 7 * (N_ + 1); ++i) {
            x_lower[i] = rmin_;
            x_upper[i] = rmax_;
        }
        x_lower[6 * (N_ + 1)] = x_upper[6 * (N_ + 1)] = r0_;

        // r
        for (int i = 7 * (N_ + 1) + 1; i < 8 * (N_ + 1); ++i) {
            x_lower[i] = xmin_;
            x_upper[i] = xmax_;
        }
        x_lower[7 * (N_ + 1)] = x_upper[7 * (N_ + 1)] = x0_;
        //x_lower[7 * (N_ + 1)] = x_upper[7 * (N_ + 1)] = r0_;

        // control input
        for (int i = 8 * (N_ + 1); i < 9 * (N_ + 1); ++i) { // 8 * (N_ + 1) + 1
            x_lower[i] = ymin_;
            x_upper[i] = ymax_;
        }
        x_lower[8 * (N_ + 1)] = x_upper[8 * (N_ + 1)] = y0_;
        //x_lower[8 * (N_ + 1)] = x_upper[8 * (N_ + 1)] = delta_v0_;

        for (int i = 9 * (N_ + 1); i < 10 * (N_ + 1); ++i) { // 9 * (N_ + 1) + 1
            x_lower[i] = delta_v_min_;
            x_upper[i] = delta_v_max_;
        }
        //x_lower[9 * (N_ + 1)] = x_upper[9 * (N_ + 1)] = delta_m0_;

        for (int i = 10 * (N_ + 1); i < 11 * (N_ + 1); ++i) { // 10 * (N_ + 1) + 1
            x_lower[i] = delta_m_min_;
            x_upper[i] = delta_m_max_;
        }
        //x_lower[10 * (N_ + 1)] = x_upper[10 * (N_ + 1)] = delta_s0_;

        for (int i = 11 * (N_ + 1); i < 12 * (N_ + 1); ++i) { // 11 * (N_ + 1) + 1
            x_lower[i] = delta_s_min_;
            x_upper[i] = delta_s_max_;
        }
        //x_lower[11 * (N_ + 1)] = x_upper[10 * (N_ + 1)] = delta_h0_;

        for (int i = 12 * (N_ + 1); i < 13 * (N_ + 1); ++i) { // 11 * (N_ + 1) + 1
            x_lower[i] = delta_h_min_;
            x_upper[i] = delta_h_max_;
        }
        //x_lower[11 * (N_ + 1)] = x_upper[10 * (N_ + 1)] = delta_h0_;

        for (int i = 13 * (N_ + 1); i < 14 * (N_ + 1); ++i) { // 11 * (N_ + 1) + 1
            x_lower[i] = delta_n_min_;
            x_upper[i] = delta_n_max_;
        }
        //x_lower[11 * (N_ + 1)] = x_upper[10 * (N_ + 1)] = delta_h0_;

        x_lower[n - 1] = tf_;
        x_upper[n - 1] = tf_;
        std::copy(x_lower.begin(), x_lower.end(), x_l);
        std::copy(x_upper.begin(), x_upper.end(), x_u);

        std::fill(g_l, g_l + m, 0);
        std::fill(g_u, g_u + m, 0);

        // for (int i = 0; i < N_; ++i) {
        //     g_l[9 * (N_ + 1) + i] = -std::numeric_limits<double>::infinity();
        //     g_u[9 * (N_ + 1) + i] = 0.0;
        // }


        return true;
    }

    virtual bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda) {
        for(Index i=0; i<n-1; ++i) 
            x[i] = 1.0;
            x[n - 1] = tf_;
        return true;
    }

    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
        const double wz = 100;//20.0;
        const double wtheta = 1;//20.0;
        const double wy = 20.0;
        const double wx = 20.0;
        const double wpsi = 10;//20.0;

        const double wdv = 100;//0.001;
        const double wdm = 0.001;//0.0001;
        const double wds = 0.12;//0.01;
        const double wdh = 0.01;//0.01;
        const double wdn = 1.0;//0.01;

        obj_value = 0.0;

        double tf_1 = x[n - 1];

        std::vector<double> zf_vector(N_ + 1, zf_);
        std::vector<double> thetaf_vector(N_ + 1, thetaf_);
        //std::vector<double> xf_vector(N_ + 1, xf_);
        //std::vector<double> yf_vector(N_ + 1, yf_);
        std::vector<double> psif_vector(N_ + 1, psif_);
        

        std::vector<double> z_vector(x, x + (N_ + 1));
        std::vector<double> theta_vector(x + 1 * (N_ + 1), x + 2 * (N_ + 1));
        //std::vector<double> x_vector(x + 7 * (N_ + 1), x + 8 * (N_ + 1));
        //std::vector<double> y_vector(x + 8 * (N_ + 1), x + 9 * (N_ + 1));
        std::vector<double> psi_vector(x + 5 * (N_ + 1), x + 6 * (N_ + 1));

        std::vector<double> dv_vector(x + 9 * (N_ + 1), x + 10 * (N_ + 1));
        std::vector<double> dm_vector(x + 10 * (N_ + 1), x + 11 * (N_ + 1));
        std::vector<double> ds_vector(x + 11 * (N_ + 1), x + 12 * (N_ + 1));
        std::vector<double> dh_vector(x + 12 * (N_ + 1), x + 13 * (N_ + 1));
        std::vector<double> dn_vector(x + 13 * (N_ + 1), x + 14 * (N_ + 1));

        std::vector<double> z_diff(N_ + 1);
        std::vector<double> theta_diff(N_ + 1);
        std::vector<double> y_diff(N_ + 1);
        std::vector<double> psi_diff(N_ + 1);

        std::vector<double> dv_sqr(N_ + 1);
        std::vector<double> dm_sqr(N_ + 1);
        std::vector<double> ds_sqr(N_ + 1);
        std::vector<double> dh_sqr(N_ + 1);
        std::vector<double> dn_sqr(N_ + 1);

        vdSub(N_ + 1, z_vector.data(), zf_vector.data(), z_diff.data());
        vdSub(N_ + 1, theta_vector.data(), thetaf_vector.data(), theta_diff.data());
        //vdSub(N_ + 1, y_vector.data(), yf_vector.data(), y_diff.data());
        vdSub(N_ + 1, psi_vector.data(), psif_vector.data(), psi_diff.data());

        vdSqr(N_ + 1, z_diff.data(), z_diff.data());
        vdSqr(N_ + 1, theta_diff.data(), theta_diff.data());
        vdSqr(N_ + 1, psi_diff.data(), psi_diff.data());
        vdSqr(N_ + 1, y_diff.data(), y_diff.data());

        vdSqr(N_ + 1, dv_vector.data(), dv_sqr.data());
        vdSqr(N_ + 1, dm_vector.data(), dm_sqr.data());
        vdSqr(N_ + 1, ds_vector.data(), ds_sqr.data());
        vdSqr(N_ + 1, dh_vector.data(), dh_sqr.data());
        vdSqr(N_ + 1, dn_vector.data(), dn_sqr.data());

        double z_diff_sum = cblas_dasum(N_ + 1, z_diff.data(), 1);
        double theta_diff_sum = cblas_dasum(N_ + 1, theta_diff.data(), 1);
        double y_diff_sum = cblas_dasum(N_ + 1, y_diff.data(), 1);
        double psi_diff_sum = cblas_dasum(N_ + 1, psi_diff.data(), 1);

        // Print x_vector:
        // std::cout << "x_vector = [";
        // for (int i = 0; i <= N_; ++i) {
        //     std::cout << x_vector[i];
        //     if (i < N_) std::cout << ", ";
        // }
        // std::cout << "]\n";

        // // Print xf_vector:
        // std::cout << "xf_vector = [";
        // for (int i = 0; i <= N_; ++i) {
        //     std::cout << xf_vector[i];
        //     if (i < N_) std::cout << ", ";
        // }
        // std::cout << "]\n";

        //double x_diff         = x_vector[N_] - xf_vector[N_];
        // std::cout << "x_vector[N_]   = " << x_vector[N_]   << "\n";
        // std::cout << "xf_vector[N_]  = " << xf_vector[N_]  << "\n";
        //double x_diff_squared = x_diff * x_diff;
        


        double dv_sum = cblas_dasum(N_ + 1, dv_sqr.data(), 1);
        double dm_sum = cblas_dasum(N_ + 1, dm_sqr.data(), 1);
        double ds_sum = cblas_dasum(N_ + 1, ds_sqr.data(), 1);
        double dh_sum = cblas_dasum(N_ + 1, dh_sqr.data(), 1);
        double dn_sum = cblas_dasum(N_ + 1, dn_sqr.data(), 1);


    
        //obj_value = wz * z_diff_sum + wtheta * theta_diff_sum + wpsi  * psi_diff_sum + wy * y_diff_sum + wx * x_diff_squared + tf_1 + wdv * dv_sum + wdm * dm_sum + wds * ds_sum + wdh * dh_sum + wdn * dn_sum; //

        obj_value = wz * z_diff_sum + wtheta * theta_diff_sum + wpsi  * psi_diff_sum + wdv * dv_sum + wdm * dm_sum + wds * ds_sum + wdh * dh_sum + wdn * dn_sum; //  

        // obj_value = wz*z_diff_sum
        //   + wtheta*theta_diff_sum
        //   + wpsi*psi_diff_sum
        //   + wy*y_diff_sum
        //   + wx*x_diff_squared
        //   + tf_1*(N_+1);

        return true;
    }

    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
        Bebot Bebot(N_, x[n - 1]);
        Bebot.calculate();
        const auto& Dm = Bebot.getDifferentiationMatrix();
        
        //std::cout << "tf (x[n-1]) = " << x[n - 1] << std::endl;

        // // Assume Dm is a flat vector of length (N_+1)*(N_+1).  Print it as an (N_+1)x(N_+1) matrix:
        // std::cout << "Dm = " << std::endl;
        // for (int i = 0; i < N_ + 1; ++i) {
        //     for (int j = 0; j < N_ + 1; ++j) {
        //         std::cout << Dm[i * (N_ + 1) + j] << "  ";
        //     }
        //     std::cout << std::endl;
        // }

        std::vector<double> z_vector(x, x + (N_ + 1));
        std::vector<double> theta_vector(x + 1 * (N_ + 1), x + 2 * (N_ + 1));
        std::vector<double> w_vector(x + 2 * (N_ + 1), x + 3 * (N_ + 1));
        std::vector<double> q_vector(x + 3 * (N_ + 1), x + 4 * (N_ + 1));
        std::vector<double> u_vector(x + 4 * (N_ + 1), x + 5 * (N_ + 1));
        std::vector<double> psi_vector(x + 5 * (N_ + 1), x + 6 * (N_ + 1));
        std::vector<double> r_vector(x + 6 * (N_ + 1), x + 7 * (N_ + 1));
        std::vector<double> x_vector(x + 7 * (N_ + 1), x + 8 * (N_ + 1));
        std::vector<double> y_vector(x + 8 * (N_ + 1), x + 9 * (N_ + 1));

        std::vector<double> delta_v_vector(x + 9 * (N_ + 1), x + 10 * (N_ + 1));
        std::vector<double> delta_m_vector(x + 10 * (N_ + 1), x + 11 * (N_ + 1));
        std::vector<double> delta_s_vector(x + 11 * (N_ + 1), x + 12 * (N_ + 1));
        std::vector<double> delta_h_vector(x + 12 * (N_ + 1), x + 13 * (N_ + 1));
        std::vector<double> delta_n_vector(x + 13 * (N_ + 1), x + 14 * (N_ + 1));

        
        std::vector<double> dyn1(N_ + 1);
        std::vector<double> dyn2(N_ + 1);
        std::vector<double> dyn3(N_ + 1);
        std::vector<double> dyn4(N_ + 1);
        std::vector<double> dyn5(N_ + 1);
        std::vector<double> dyn6(N_ + 1);
        std::vector<double> dyn7(N_ + 1);
        std::vector<double> dyn8(N_ + 1);
        std::vector<double> dyn9(N_ + 1);

        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, z_vector.data(), 1, 0.0, dyn1.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, theta_vector.data(), 1, 0.0, dyn2.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, w_vector.data(), 1, 0.0, dyn3.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, q_vector.data(), 1, 0.0, dyn4.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, u_vector.data(), 1, 0.0, dyn5.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, psi_vector.data(), 1, 0.0, dyn6.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, r_vector.data(), 1, 0.0, dyn7.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, x_vector.data(), 1, 0.0, dyn8.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, y_vector.data(), 1, 0.0, dyn9.data(), 1);

        std::vector<double> X1_matrix_flat(4 * (N_ + 1));
        std::vector<double> U_matrix_flat(3 * (N_ + 1));

        for (Index i = 0; i < (N_ + 1); ++i) {
            X1_matrix_flat[i] = z_vector[i];
            X1_matrix_flat[(N_ + 1) + i] = theta_vector[i];
            X1_matrix_flat[2 * (N_ + 1) + i] = w_vector[i];
            X1_matrix_flat[3 * (N_ + 1) + i] = q_vector[i];
        }

        for (Index i = 0; i < (N_ + 1); ++i) {
            U_matrix_flat[i] = delta_v_vector[i];
            U_matrix_flat[(N_ + 1) + i] = delta_m_vector[i];
            U_matrix_flat[2 * (N_ + 1) + i] = delta_s_vector[i];
        }

        std::vector<double> X2_matrix_flat(4 * (N_ + 1), 0.0);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4, (N_ + 1), 4, 1.0, &A_[0][0], 4, X1_matrix_flat.data(), (N_ + 1), 0.0, X2_matrix_flat.data(), (N_ + 1));
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4, (N_ + 1), 3, 1.0, &B_[0][0], 3, U_matrix_flat.data(), (N_ + 1), 1.0, X2_matrix_flat.data(), (N_ + 1));

        std::vector<double> X1_matrix_flat_cd(3 * (N_ + 1));
        std::vector<double> U_matrix_flat_cd(2 * (N_ + 1));

        for (Index i = 0; i < (N_ + 1); ++i) {
            X1_matrix_flat_cd[i] = u_vector[i];
            X1_matrix_flat_cd[(N_ + 1) + i] = psi_vector[i];
            X1_matrix_flat_cd[2 * (N_ + 1) + i] = r_vector[i];
        }

        for (Index i = 0; i < (N_ + 1); ++i) {
            U_matrix_flat_cd[i] = delta_h_vector[i];
            U_matrix_flat_cd[(N_ + 1) + i] = delta_n_vector[i];
        }

        std::vector<double> X2_matrix_flat_cd(3 * (N_ + 1), 0.0);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, (N_ + 1), 3, 1.0, &C_[0][0], 3, X1_matrix_flat_cd.data(), (N_ + 1), 0.0, X2_matrix_flat_cd.data(), (N_ + 1));
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, (N_ + 1), 2, 1.0, &D_[0][0], 2, U_matrix_flat_cd.data(), (N_ + 1), 1.0, X2_matrix_flat_cd.data(), (N_ + 1));



        std::vector<double> z_x2(N_ + 1);
        std::vector<double> theta_x2(N_ + 1);
        std::vector<double> w_x2(N_ + 1);
        std::vector<double> q_x2(N_ + 1);
        std::vector<double> u_x2(N_ + 1);
        std::vector<double> psi_x2(N_ + 1);
        std::vector<double> r_x2(N_ + 1);

        std::vector<double> x_x2(N_ + 1);
        std::vector<double> y_x2(N_ + 1);

        for (Index i = 0; i < (N_ + 1); ++i) {
            z_x2[i] = X2_matrix_flat[i];
            theta_x2[i] = X2_matrix_flat[(N_ + 1) + i];
            w_x2[i] = X2_matrix_flat[2 * (N_ + 1) + i];
            q_x2[i] = X2_matrix_flat[3 * (N_ + 1) + i];
        }

        for (Index i = 0; i < (N_ + 1); ++i) {
            u_x2[i] = X2_matrix_flat_cd[i];
            psi_x2[i] = X2_matrix_flat_cd[(N_ + 1) + i];
            r_x2[i] = X2_matrix_flat_cd[2 * (N_ + 1) + i];
        }

        // // --- 1) Print psi_vector before any multiplication ---
        // std::cout << "psi_vector (length " << (N_ + 1) << "): [";
        // for (int i = 0; i < N_ + 1; ++i) {
        //     std::cout << psi_vector[i];
        //     if (i < N_) std::cout << ", ";
        // }
        // std::cout << "]\n\n";
        // std::cout << "i |   psi_i    |  cos(psi_i)  |  sin(psi_i)  |   x_x2[i]    |   y_x2[i]\n";
        // std::cout << "---------------------------------------------------------------------\n";
        for (Index i = 0; i < N_ + 1; ++i) {
            double u_i    = delta_n_vector[i] * 4.0;
            double cospsi_i = std::cos(psi_vector[i]);
            double sinpsi_i = std::sin(psi_vector[i]);
            x_x2[i] = u_vector[i] * cospsi_i;
            y_x2[i] = u_vector[i] * sinpsi_i;

            // // print row‐by‐row
            // std::cout << std::setw(2) << i << " | "
            //         << std::setw(9)  << psi_vector[i] << " | "
            //         << std::setw(11) << cospsi_i       << " | "
            //         << std::setw(11) << sinpsi_i       << " | "
            //         << std::setw(11) << x_x2[i]         << " | "
            //         << std::setw(11) << y_x2[i]         << "\n";
        }

        // std::cout << "\n";
        

        std::vector<double> g1(N_ + 1);
        std::vector<double> g2(N_ + 1);
        std::vector<double> g3(N_ + 1);
        std::vector<double> g4(N_ + 1);
        std::vector<double> g5(N_ + 1);
        std::vector<double> g6(N_ + 1);
        std::vector<double> g7(N_ + 1);
        std::vector<double> g8(N_ + 1);
        std::vector<double> g9(N_ + 1);

        vdSub(N_ + 1, dyn1.data(), z_x2.data(), g1.data());
        vdSub(N_ + 1, dyn2.data(), theta_x2.data(), g2.data());
        vdSub(N_ + 1, dyn3.data(), w_x2.data(), g3.data());
        vdSub(N_ + 1, dyn4.data(), q_x2.data(), g4.data());
        vdSub(N_ + 1, dyn5.data(), u_x2.data(), g5.data());
        vdSub(N_ + 1, dyn6.data(), psi_x2.data(), g6.data());
        vdSub(N_ + 1, dyn7.data(), r_x2.data(), g7.data());
        vdSub(N_ + 1, dyn8.data(), x_x2.data(), g8.data());
        vdSub(N_ + 1, dyn9.data(), y_x2.data(), g9.data());

        for (Index i = 0; i < (N_ + 1); ++i) {
            g[i] = g1[i];
            g[(N_ + 1) + i] = g2[i];
            g[2 * (N_ + 1) + i] = g3[i];
            g[3 * (N_ + 1) + i] = g4[i];
            g[4 * (N_ + 1) + i] = g5[i];
            g[5 * (N_ + 1) + i] = g6[i];
            g[6 * (N_ + 1) + i] = g7[i];
            g[7 * (N_ + 1) + i] = g8[i];
            g[8 * (N_ + 1) + i] = g9[i];
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


        // Ensure solution_x is 7 * (N + 1)
        solution_x_.resize(14 * (N_ + 1) + 1);
        
        // Copy the first 4 * (N + 1) elements from the x array
        for (Index i = 0; i < 14 * (N_ + 1) + 1; ++i) {
            solution_x_[i] = x[i];
        }
        
        

        // Copy the remaining 3 * (N + 1) elements from the x array
        //for (Index i = 0; i < 2 * (N_ + 1); ++i) {
        //    solution_x_[4* (N_ + 1) + i] = x[4 * (N_ + 1) + i]; //-
        //}

        // Copy the remaining 3 * (N + 1) elements from the x array
        for (Index i = 0; i < 1 * (N_ + 1); ++i) {
            solution_x_[11* (N_ + 1) + i] = x[11 * (N_ + 1) + i]; //-
        }
        
        double tf_ = solution_x_.back();

        final_obj_value_ = obj_value;         
        bebot_ = Bebot(N_, tf_);
        bebot_.calculate();
        
        final_time_.resize(1000);
        for (int i = 0; i < 1000; ++i) {
            final_time_[i] = i * tf_ / 999.0;
        }

        std::vector<double> z_vector(solution_x_.begin(), solution_x_.begin() + (N_ + 1));
        std::vector<double> theta_vector(solution_x_.begin() + (N_ + 1), solution_x_.begin() + 2 * (N_ + 1));
        std::vector<double> w_vector(solution_x_.begin() + 2 * (N_ + 1), solution_x_.begin() + 3 * (N_ + 1));
        std::vector<double> q_vector(solution_x_.begin() + 3 * (N_ + 1), solution_x_.begin() + 4 * (N_ + 1));
        std::vector<double> u_vector(solution_x_.begin() + 4 * (N_ + 1), solution_x_.begin() + 5 * (N_ + 1));      
        std::vector<double> psi_vector(solution_x_.begin() + 5 * (N_ + 1), solution_x_.begin() + 6 * (N_ + 1));
        std::vector<double> r_vector(solution_x_.begin() + 6 * (N_ + 1), solution_x_.begin() + 7 * (N_ + 1));

        std::vector<double> x_vector(solution_x_.begin() + 7 * (N_ + 1), solution_x_.begin() + 8 * (N_ + 1));
        std::vector<double> y_vector(solution_x_.begin() + 8 * (N_ + 1), solution_x_.begin() + 9 * (N_ + 1));

        std::vector<double> delta_v_vector(solution_x_.begin() + 9 * (N_ + 1), solution_x_.begin() + 10 * (N_ + 1));
        std::vector<double> delta_m_vector(solution_x_.begin() + 10 * (N_ + 1), solution_x_.begin() + 11 * (N_ + 1));
        std::vector<double> delta_s_vector(solution_x_.begin() + 11 * (N_ + 1), solution_x_.begin() + 12 * (N_ + 1));
        std::vector<double> delta_h_vector(solution_x_.begin() + 12 * (N_ + 1), solution_x_.begin() + 13 * (N_ + 1));
        std::vector<double> delta_n_vector(solution_x_.begin() + 13 * (N_ + 1), solution_x_.end() - 1);
        

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
        // print_vec("y_vector",        y_vector);
        // print_vec("psi_vector",      psi_vector);
        // print_vec("v_vector",        v_vector);
        // print_vec("r_vector",        r_vector);

        // print_vec("delta_v_vector",  delta_v_vector);
        // print_vec("delta_m_vector",  delta_m_vector);
        // print_vec("delta_s_vector",  delta_s_vector);
        // print_vec("delta_h_vector",  delta_h_vector);


        std::vector<std::vector<double>> z_2d(1, z_vector);
        std::vector<std::vector<double>> theta_2d(1, theta_vector);
        std::vector<std::vector<double>> w_2d(1, w_vector);
        std::vector<std::vector<double>> q_2d(1, q_vector);

        std::vector<std::vector<double>> u_2d(1, u_vector);
        std::vector<std::vector<double>> psi_2d(1, psi_vector);
        std::vector<std::vector<double>> r_2d(1, r_vector);

        std::vector<std::vector<double>> x_2d(1, x_vector);
        std::vector<std::vector<double>> y_2d(1, y_vector);

        std::vector<std::vector<double>> delta_v_2d(1, delta_v_vector);
        std::vector<std::vector<double>> delta_m_2d(1, delta_m_vector);
        std::vector<std::vector<double>> delta_s_2d(1, delta_s_vector);
        std::vector<std::vector<double>> delta_h_2d(1, delta_h_vector);
        std::vector<std::vector<double>> delta_n_2d(1, delta_n_vector);

        std::vector<std::vector<double>> bernstein_z = BernsteinPoly(z_2d, final_time_, 0, tf_);
        std::vector<std::vector<double>> bernstein_theta = BernsteinPoly(theta_2d, final_time_, 0, tf_);
        std::vector<std::vector<double>> bernstein_w = BernsteinPoly(w_2d, final_time_, 0, tf_);
        std::vector<std::vector<double>> bernstein_q = BernsteinPoly(q_2d, final_time_, 0, tf_);

        std::vector<std::vector<double>> bernstein_u = BernsteinPoly(u_2d, final_time_, 0, tf_);
        std::vector<std::vector<double>> bernstein_psi = BernsteinPoly(psi_2d, final_time_, 0, tf_);
        std::vector<std::vector<double>> bernstein_r = BernsteinPoly(r_2d, final_time_, 0, tf_);

        std::vector<std::vector<double>> bernstein_x = BernsteinPoly(x_2d, final_time_, 0, tf_);
        std::vector<std::vector<double>> bernstein_y = BernsteinPoly(y_2d, final_time_, 0, tf_);

        std::vector<std::vector<double>> bernstein_delta_v = BernsteinPoly(delta_v_2d, final_time_, 0, tf_);
        std::vector<std::vector<double>> bernstein_delta_m = BernsteinPoly(delta_m_2d, final_time_, 0, tf_);
        std::vector<std::vector<double>> bernstein_delta_s = BernsteinPoly(delta_s_2d, final_time_, 0, tf_);
        std::vector<std::vector<double>> bernstein_delta_h = BernsteinPoly(delta_h_2d, final_time_, 0, tf_);
        std::vector<std::vector<double>> bernstein_delta_n = BernsteinPoly(delta_n_2d, final_time_, 0, tf_);

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


        writeToCSV(final_time_, flatten(bernstein_u), "u.csv");
        writeToCSV(bebot_.getNodes(), u_vector, "u_controlpoints.csv");
        writeToCSV(final_time_, flatten(bernstein_psi), "psi.csv");
        writeToCSV(bebot_.getNodes(), psi_vector, "psi_controlpoints.csv");
        writeToCSV(final_time_, flatten(bernstein_r), "r.csv");
        writeToCSV(bebot_.getNodes(), r_vector, "r_controlpoints.csv");


        writeToCSV(final_time_, flatten(bernstein_x), "x.csv");
        writeToCSV(bebot_.getNodes(), x_vector, "x_controlpoints.csv");
        writeToCSV(final_time_, flatten(bernstein_y), "y.csv");
        writeToCSV(bebot_.getNodes(), y_vector, "y_controlpoints.csv");

        writeToCSV(final_time_, flatten(bernstein_delta_v), "delta_v.csv");
        writeToCSV(bebot_.getNodes(), delta_v_vector, "delta_v_controlpoints.csv");
        writeToCSV(final_time_, flatten(bernstein_delta_m), "delta_m.csv");
        writeToCSV(bebot_.getNodes(), delta_m_vector, "delta_m_controlpoints.csv");
        writeToCSV(final_time_, flatten(bernstein_delta_s), "delta_s.csv");
        writeToCSV(bebot_.getNodes(), delta_s_vector, "delta_s_controlpoints.csv");
        writeToCSV(final_time_, flatten(bernstein_delta_h), "delta_h.csv");
        writeToCSV(bebot_.getNodes(), delta_h_vector, "delta_h_controlpoints.csv");
        writeToCSV(final_time_, flatten(bernstein_delta_n), "delta_n.csv");
        writeToCSV(bebot_.getNodes(), delta_n_vector, "delta_n_controlpoints.csv");
        std::cout << "Solution finalized and written to CSV files" << std::endl;

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
    double delta_n_max_;
    double delta_n_min_;
    double zmax_;
    double zmin_;
    double wmax_;  
    double wmin_;
    double thetamax_;
    double thetamin_;
    double qmax_;
    double qmin_;
    double umax_;
    double umin_;
    double psimax_;  
    double psimin_;
    double rmax_;
    double rmin_;
    double xmax_;
    double xmin_;
    double ymax_;
    double ymin_;
    double z0_;
    double w0_;
    double theta0_;
    double q0_;
    
    double u0_;
    double psi0_;
    double r0_;
    double x0_;
    double y0_;
    double delta_v0_;
    double delta_s0_;
    double delta_m0_;
    double delta_h0_;
    double delta_n0_;
    double zf_;
    double thetaf_;
    double xf_;
    double yf_;
    double psif_;
    double t0_;
    double tend_;

    std::array<std::array<double, 4>, 4> A_;
    std::array<std::array<double, 3>, 4> B_;
    std::array<std::array<double, 3>, 3> C_;
    std::array<std::array<double, 2>, 3> D_;
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
        double delta_s_max, double delta_s_min, double delta_m_max, double delta_m_min, double delta_h_max, double delta_h_min, double delta_n_max, double delta_n_min,
        double zmax, double zmin, double wmax, double wmin, double thetamax, double thetamin, double qmax, double qmin, 
        double umax, double umin, double psimax, double psimin, double rmax, double rmin, 
        double xmax, double xmin, double ymax, double ymin, 
        double z0, double w0, double theta0, double q0,
        double u0, double psi0, double r0, 
        double x0, double y0, 
        double delta_v0, double delta_s0, double delta_m0, double delta_h0, double delta_n0,
        double zf, double thetaf, double xf, double yf, double psif, 
        double a11, double a12, double a13, double a14, 
        double a21, double a22, double a23, double a24, 
        double a31, double a32, double a33, double a34, 
        double a41, double a42, double a43, double a44, 
        double b11, double b12, double b13, 
        double b21, double b22, double b23, 
        double b31, double b32, double b33, 
        double b41, double b42, double b43, 
        double c11, double c12, double c13, 
        double c21, double c22, double c23, 
        double c31, double c32, double c33, 
        double d11, double d12,  
        double d21, double d22, 
        double d31, double d32,
        double t0, double tend) {
        
        std::cout << std::fixed << std::setprecision(5);
        //Print out every incoming argument:
        std::cout << "==== create_point_set_problem called ====\n";
        std::cout << "N = " << N << "\n";
        std::cout << "tf = " << tf << "\n";

        std::cout << "delta_v_max = " << delta_v_max
                  << ", delta_v_min = " << delta_v_min << "\n";
        std::cout << "delta_s_max = " << delta_s_max
                  << ", delta_s_min = " << delta_s_min << "\n";
        std::cout << "delta_m_max = " << delta_m_max
                  << ", delta_m_min = " << delta_m_min << "\n";
        std::cout << "delta_h_max = " << delta_h_max
                  << ", delta_h_min = " << delta_h_min << "\n";
        std::cout << "delta_n_max = " << delta_n_max
                  << ", delta_n_min = " << delta_n_min << "\n\n";

        std::cout << "zmax = " << zmax << ", zmin = " << zmin << "\n";
        std::cout << "wmax = " << wmax << ", wmin = " << wmin << "\n";
        std::cout << "thetamax = " << thetamax << ", thetamin = " << thetamin << "\n";
        std::cout << "qmax = " << qmax << ", qmin = " << qmin << "\n\n";

        std::cout << "umax = " << umax << ", umin = " << umin << "\n";
        std::cout << "psimax = " << psimax << ", psimin = " << psimin << "\n";
        std::cout << "rmax = " << rmax << ", rmin = " << rmin << "\n\n";

        std::cout << "xmax = " << xmax << ", xmin = " << xmin << "\n";
        std::cout << "ymax = " << ymax << ", ymin = " << ymin << "\n\n";

        std::cout << "Initial states:\n";
        std::cout << "  z0 = " << z0 << ", w0 = " << w0
                  << ", theta0 = " << theta0 << ", q0 = " << q0 << "\n";
        std::cout << "  u0 = " << u0 << ", psi0 = " << psi0 << ", r0 = " << r0 << "\n";
        std::cout << "  x0 = " << x0 << ", y0 = " << y0 << "\n\n";

        std::cout << "Initial controls:\n";
        std::cout << "  delta_v0 = " << delta_v0
                  << ", delta_s0 = " << delta_s0
                  << ", delta_m0 = " << delta_m0
                  << ", delta_h0 = " << delta_h0
                  << ", delta_n0 = " << delta_n0 << "\n\n";

        std::cout << "Final targets:\n";
        std::cout << "  zf = " << zf
                  << ", thetaf = " << thetaf
                  << ", xf = " << xf
                  << ", yf = " << yf
                  << ", psif = " << psif << "\n\n";

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

        // std::cout << "C‐matrix (3×3):\n";
        // std::cout << "  [" << c11 << ", " << c12 << ", " << c13 << "]\n";
        // std::cout << "  [" << c21 << ", " << c22 << ", " << c23 << "]\n";
        // std::cout << "  [" << c31 << ", " << c32 << ", " << c33 << "]\n\n";

        // std::cout << "D‐matrix (3×2):\n";
        // std::cout << "  [" << d11 << ", " << d12 << "]\n";
        // std::cout << "  [" << d21 << ", " << d22 << "]\n";
        // std::cout << "  [" << d31 << ", " << d32 << "]\n\n";

        std::cout << "t0 = " << t0 << ", tend = " << tend << "\n";
        std::cout << "===========================================\n\n";
        return new PointSetProblem(N, tf, delta_v_max, delta_v_min, delta_s_max, delta_s_min, 
            delta_m_max, delta_m_min, delta_h_max, delta_h_min, delta_n_max, delta_n_min,
            zmax, zmin, wmax, wmin, thetamax, thetamin, qmax, qmin, 
            umax, umin, psimax, psimin, rmax, rmin,
            xmax, xmin, ymax, ymin, 
            z0, w0, theta0, q0, u0, psi0, r0,
            x0, y0, 
            delta_v0, delta_s0, delta_m0, delta_h0, delta_n0,
            zf, thetaf, xf, yf, psif,
            a11, a12, a13, a14, 
            a21, a22, a23, a24, 
            a31, a32, a33, a34, 
            a41, a42, a43, a44, 
            b11, b12, b13, 
            b21, b22, b23, 
            b31, b32, b33, 
            b41, b42, b43, 
            c11, c12, c13, 
            c21, c22, c23, 
            c31, c32, c33, 
            d11, d12,
            d21, d22, 
            d31, d32,            
            t0, tend);
    }

    void solve_point_set_problem(PointSetProblem* problem) {
        SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
        app->Options()->SetStringValue("linear_solver", "ma57");
        app->Options()->SetStringValue("mu_strategy", "adaptive");
        app->Options()->SetStringValue("gradient_approximation", "finite-difference-values");
        app->Options()->SetStringValue("jacobian_approximation", "finite-difference-values");
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");
        app->Options()->SetIntegerValue("max_iter", 400);
        app->Options()->SetNumericValue("tol",             1e-3);   // OptimalityTolerance = 1e-3
        app->Options()->SetNumericValue("constr_viol_tol", 1e-6);
        app->Options()->SetNumericValue("acceptable_tol",        1e-6);
        // somewhere before app->Initialize():
        //app->Options()->SetIntegerValue("max_line_search_step_retries", 200);
        //app->Options()->SetNumericValue("alpha_for_y", 0.6);
        //app->Options()->SetNumericValue("beta_for_y",  0.4);


        //app->Options()->SetNumericValue("constr_viol_tol", 1e-6);
        app->Options()->SetIntegerValue("print_level", 0); 
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
