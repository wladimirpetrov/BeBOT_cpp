// #include "../../../../Ipopt_ma57_solver/src/Interfaces/IpIpoptApplication.hpp"
// #include "../../../../Ipopt_ma57_solver/src/Interfaces/IpTNLP.hpp"

// #include <cmath>
// #include <iostream>
// #include <fstream>
// #include <vector>
// #include <iomanip>
// #include <array>
// #include <string>
// #include <limits>
// #include <algorithm>
// #include <sstream>

// #include "../../../../include/bebot.h"
// #include "../../../../include/bernsteinpoly.h"
// #include "../../../../include/bernsteinproduct.h"
// #include "../../../../include/degelevmatrix.h"

// #include "mkl.h"

// using namespace Ipopt;

// class PointSetProblem : public Ipopt::TNLP {
// public:
//     PointSetProblem(
//         int N,
//         double tf,
//         double px_max,
//         double px_min,
//         double py_max,
//         double py_min,
//         double pz_max,
//         double pz_min,
//         double psi_max,
//         double psi_min,
//         double v_max,
//         double w_max,
//         double a_max,
//         double aw_max,
//         double px_cur,
//         double py_cur,
//         double pz_cur,
//         double psi_cur,
//         double vx_cur,
//         double vy_cur,
//         double vz_cur,
//         double w_cur,
//         double pxf,
//         double pyf,
//         double pzf,
//         double psif
//     )
//         : N_(N),
//           tf_(tf),
//           px_max_(px_max),
//           px_min_(px_min),
//           py_max_(py_max),
//           py_min_(py_min),
//           pz_max_(pz_max),
//           pz_min_(pz_min),
//           psi_max_(psi_max),
//           psi_min_(psi_min),
//           v_max_(v_max),
//           w_max_(w_max),
//           a_max_(a_max),
//           aw_max_(aw_max),
//           px_cur_(px_cur),
//           py_cur_(py_cur),
//           pz_cur_(pz_cur),
//           psi_cur_(psi_cur),
//           vx_cur_(vx_cur),
//           vy_cur_(vy_cur),
//           vz_cur_(vz_cur),
//           w_cur_(w_cur),
//           pxf_(pxf),
//           pyf_(pyf),
//           pzf_(pzf),
//           psif_(psif),
//           obs_x_{{-0.6, 0.6}},
//           obs_y_{{ 2.0, 1.8}},
//           obs_sep_(0.3),
//           obs_deg_elev_extra_(5),
//           bebot_(N, tf_) {
//         bebot_.calculate();

//         const int L = N_ + 1;

//         tau_vec_.resize(L);
//         r_tau_.resize(L);
//         r_tau_row_.resize(L);

//         if (L == 1) {
//             const double tau0 = -1.0;
//             const double denom = 1.0 - tau0;

//             tau_vec_[0] = tau0;
//             r_tau_[0] = 2.0 / (denom * denom);
//             r_tau_row_[0] = r_tau_[0];
//         } else {
//             const double tau0 = -1.0;
//             const double tau1 = 0.9;
//             const double step = (tau1 - tau0) / static_cast<double>(N_);

//             for (int i = 0; i < L; ++i) {
//                 const double tau = tau0 + step * static_cast<double>(i);
//                 const double denom = 1.0 - tau;

//                 tau_vec_[i] = tau;
//                 r_tau_[i] = 2.0 / (denom * denom);
//                 r_tau_row_[i] = r_tau_[i];
//             }
//         }
//     }

//     void writeToCSV(
//         const std::vector<double>& times,
//         const std::vector<double>& values,
//         const std::string& filename
//     ) {
//         std::ofstream outFile(filename);

//         if (!outFile.is_open()) {
//             std::cerr << "Failed to open file: " << filename << std::endl;
//             return;
//         }

//         outFile << "Time,Value\n";

//         for (size_t i = 0; i < times.size(); ++i) {
//             outFile << std::fixed << std::setprecision(6)
//                     << times[i] << "," << values[i] << "\n";
//         }

//         outFile.close();
//     }

//     virtual bool get_nlp_info(
//         Index& n,
//         Index& m,
//         Index& nnz_jac_g,
//         Index& nnz_h_lag,
//         IndexStyleEnum& index_style
//     ) {
//         const int L = N_ + 1;

//         n = 12 * L;

//         /*
//             Existing constraints:
//             8 * (N + 1): dynamics equality constraints
//             1 * (N + 1): speed inequality constraints
//             1 * (N + 1): acceleration inequality constraints

//             Obstacle constraints:
//             For each obstacle:
//                 dist2obs_square has degree 2N
//                 then degree-elevated to 2N + obs_deg_elev_extra_
//                 number of coefficients = 2N + obs_deg_elev_extra_ + 1

//             Total:
//                 m = 10 * (N + 1)
//                     + number_of_obstacles * (2N + obs_deg_elev_extra_ + 1)
//         */

//         const int obs_degree = 2 * N_ + obs_deg_elev_extra_;
//         const int obs_constraints =
//             static_cast<int>(obs_x_.size()) * (obs_degree + 1);

//         m = 10 * L + obs_constraints;

//         nnz_jac_g = n * m;
//         nnz_h_lag = 0;
//         index_style = TNLP::C_STYLE;

//         return true;
//     }

//     virtual bool get_bounds_info(
//         Index n,
//         Number* x_l,
//         Number* x_u,
//         Index m,
//         Number* g_l,
//         Number* g_u
//     ) {
//         const int L = N_ + 1;

//         std::vector<double> x_lower(
//             n,
//             -std::numeric_limits<double>::infinity()
//         );

//         std::vector<double> x_upper(
//             n,
//             std::numeric_limits<double>::infinity()
//         );

//         // px
//         for (int i = 1; i < L; ++i) {
//             x_lower[i] = px_min_;
//             x_upper[i] = px_max_;
//         }
//         x_lower[0] = px_cur_;
//         x_upper[0] = px_cur_;

//         // py
//         for (int i = L + 1; i < 2 * L; ++i) {
//             x_lower[i] = py_min_;
//             x_upper[i] = py_max_;
//         }
//         x_lower[L] = py_cur_;
//         x_upper[L] = py_cur_;

//         // pz
//         for (int i = 2 * L + 1; i < 3 * L; ++i) {
//             x_lower[i] = pz_min_;
//             x_upper[i] = pz_max_;
//         }
//         x_lower[2 * L] = pz_cur_;
//         x_upper[2 * L] = pz_cur_;

//         // psi
//         for (int i = 3 * L + 1; i < 4 * L; ++i) {
//             x_lower[i] = psi_min_;
//             x_upper[i] = psi_max_;
//         }
//         x_lower[3 * L] = psi_cur_;
//         x_upper[3 * L] = psi_cur_;

//         // vx
//         for (int i = 4 * L + 1; i < 5 * L; ++i) {
//             x_lower[i] = -v_max_;
//             x_upper[i] =  v_max_;
//         }
//         x_lower[4 * L] = vx_cur_;
//         x_upper[4 * L] = vx_cur_;

//         // vy
//         for (int i = 5 * L + 1; i < 6 * L; ++i) {
//             x_lower[i] = -v_max_;
//             x_upper[i] =  v_max_;
//         }
//         x_lower[5 * L] = vy_cur_;
//         x_upper[5 * L] = vy_cur_;

//         // vz
//         for (int i = 6 * L + 1; i < 7 * L; ++i) {
//             x_lower[i] = -v_max_;
//             x_upper[i] =  v_max_;
//         }
//         x_lower[6 * L] = vz_cur_;
//         x_upper[6 * L] = vz_cur_;

//         // w
//         for (int i = 7 * L + 1; i < 8 * L; ++i) {
//             x_lower[i] = -w_max_;
//             x_upper[i] =  w_max_;
//         }
//         x_lower[7 * L] = w_cur_;
//         x_upper[7 * L] = w_cur_;

//         // ax
//         for (int i = 8 * L; i < 9 * L; ++i) {
//             x_lower[i] = -a_max_;
//             x_upper[i] =  a_max_;
//         }

//         // ay
//         for (int i = 9 * L; i < 10 * L; ++i) {
//             x_lower[i] = -a_max_;
//             x_upper[i] =  a_max_;
//         }

//         // az
//         for (int i = 10 * L; i < 11 * L; ++i) {
//             x_lower[i] = -a_max_;
//             x_upper[i] =  a_max_;
//         }

//         // aw
//         for (int i = 11 * L; i < 12 * L; ++i) {
//             x_lower[i] = -aw_max_;
//             x_upper[i] =  aw_max_;
//         }

//         std::copy(x_lower.begin(), x_lower.end(), x_l);
//         std::copy(x_upper.begin(), x_upper.end(), x_u);

//         // Dynamics equality constraints
//         for (int i = 0; i < 8 * L; ++i) {
//             g_l[i] = 0.0;
//             g_u[i] = 0.0;
//         }

//         // Speed inequality constraints: c_speed <= 0
//         for (int i = 8 * L; i < 9 * L; ++i) {
//             g_l[i] = -std::numeric_limits<double>::infinity();
//             g_u[i] = 0.0;
//         }

//         // Acceleration inequality constraints: c_accel <= 0
//         for (int i = 9 * L; i < 10 * L; ++i) {
//             g_l[i] = -std::numeric_limits<double>::infinity();
//             g_u[i] = 0.0;
//         }

//         // Obstacle inequality constraints: c_obs <= 0
//         const int obs_start = 10 * L;

//         for (Index i = obs_start; i < m; ++i) {
//             g_l[i] = -std::numeric_limits<double>::infinity();
//             g_u[i] = 0.0;
//         }

//         return true;
//     }

//     virtual bool get_starting_point(
//         Index n,
//         bool init_x,
//         Number* x,
//         bool init_z,
//         Number* z_L,
//         Number* z_U,
//         Index m,
//         bool init_lambda,
//         Number* lambda
//     ) {
//         for (Index i = 0; i < n; ++i) {
//             x[i] = 1.0;
//         }

//         return true;
//     }

//     virtual bool eval_f(
//         Index n,
//         const Number* x,
//         bool new_x,
//         Number& obj_value
//     ) {
//         const int L = N_ + 1;

//         const double w_p   = 0.01;
//         const double w_psi = 1.0;
//         const double w_a   = 0.01;
//         const double w_aw  = 0.001;

//         obj_value = 0.0;

//         std::vector<double> pxf_vector(L, pxf_);
//         std::vector<double> pyf_vector(L, pyf_);
//         std::vector<double> pzf_vector(L, pzf_);
//         std::vector<double> psif_vector(L, psif_);

//         std::vector<double> px_vector(x, x + L);
//         std::vector<double> py_vector(x + 1 * L, x + 2 * L);
//         std::vector<double> pz_vector(x + 2 * L, x + 3 * L);
//         std::vector<double> psi_vector(x + 3 * L, x + 4 * L);

//         std::vector<double> ax_vector(x + 8  * L, x + 9  * L);
//         std::vector<double> ay_vector(x + 9  * L, x + 10 * L);
//         std::vector<double> az_vector(x + 10 * L, x + 11 * L);
//         std::vector<double> aw_vector(x + 11 * L, x + 12 * L);

//         std::vector<double> px_diff(L);
//         std::vector<double> py_diff(L);
//         std::vector<double> pz_diff(L);
//         std::vector<double> psi_diff(L);

//         vdSub(L, px_vector.data(),  pxf_vector.data(),  px_diff.data());
//         vdSub(L, py_vector.data(),  pyf_vector.data(),  py_diff.data());
//         vdSub(L, pz_vector.data(),  pzf_vector.data(),  pz_diff.data());
//         vdSub(L, psi_vector.data(), psif_vector.data(), psi_diff.data());

//         std::vector<double> px_diff_sqr(L);
//         std::vector<double> py_diff_sqr(L);
//         std::vector<double> pz_diff_sqr(L);
//         std::vector<double> psi_diff_sqr(L);

//         vdSqr(L, px_diff.data(),  px_diff_sqr.data());
//         vdSqr(L, py_diff.data(),  py_diff_sqr.data());
//         vdSqr(L, pz_diff.data(),  pz_diff_sqr.data());
//         vdSqr(L, psi_diff.data(), psi_diff_sqr.data());

//         std::vector<double> px_weighted(L);
//         std::vector<double> py_weighted(L);
//         std::vector<double> pz_weighted(L);
//         std::vector<double> psi_weighted(L);

//         vdMul(L, r_tau_.data(), px_diff_sqr.data(),  px_weighted.data());
//         vdMul(L, r_tau_.data(), py_diff_sqr.data(),  py_weighted.data());
//         vdMul(L, r_tau_.data(), pz_diff_sqr.data(),  pz_weighted.data());
//         vdMul(L, r_tau_.data(), psi_diff_sqr.data(), psi_weighted.data());

//         std::vector<double> ones(L, 1.0);

//         const double sum_px_state =
//             cblas_ddot(L, px_weighted.data(), 1, ones.data(), 1);

//         const double sum_py_state =
//             cblas_ddot(L, py_weighted.data(), 1, ones.data(), 1);

//         const double sum_pz_state =
//             cblas_ddot(L, pz_weighted.data(), 1, ones.data(), 1);

//         const double sum_psi_state =
//             cblas_ddot(L, psi_weighted.data(), 1, ones.data(), 1);

//         const double state_term =
//             w_p   * sum_px_state
//           + w_p   * sum_py_state
//           + w_p   * sum_pz_state
//           + w_psi * sum_psi_state;

//         std::vector<double> ax_sqr(L);
//         std::vector<double> ay_sqr(L);
//         std::vector<double> az_sqr(L);
//         std::vector<double> aw_sqr(L);

//         vdSqr(L, ax_vector.data(), ax_sqr.data());
//         vdSqr(L, ay_vector.data(), ay_sqr.data());
//         vdSqr(L, az_vector.data(), az_sqr.data());
//         vdSqr(L, aw_vector.data(), aw_sqr.data());

//         const double sum_ax =
//             cblas_ddot(L, ax_sqr.data(), 1, ones.data(), 1);

//         const double sum_ay =
//             cblas_ddot(L, ay_sqr.data(), 1, ones.data(), 1);

//         const double sum_az =
//             cblas_ddot(L, az_sqr.data(), 1, ones.data(), 1);

//         const double sum_aw =
//             cblas_ddot(L, aw_sqr.data(), 1, ones.data(), 1);

//         const double control_term =
//             w_a  * sum_ax
//           + w_a  * sum_ay
//           + w_a  * sum_az
//           + w_aw * sum_aw;

//         obj_value = state_term + control_term;

//         return true;
//     }

//     virtual bool eval_g(
//         Index n,
//         const Number* x,
//         bool new_x,
//         Index m,
//         Number* g
//     ) {
//         const int L = N_ + 1;

//         /*
//             CHANGED ONLY THIS PART:

//             Before, this function constructed and calculated BeBOT every time:

//                 Bebot Bebot(N_, tf_);
//                 Bebot.calculate();
//                 const auto& Dm = Bebot.getDifferentiationMatrix();

//             That is expensive because eval_g() is called many times by IPOPT.

//             Now, we reuse the already-computed member bebot_.
//             bebot_.calculate() is called once in the constructor.
//         */
//         const auto& Dm = bebot_.getDifferentiationMatrix();

//         std::vector<double> px_vector(x, x + L);
//         std::vector<double> py_vector(x + 1 * L, x + 2 * L);
//         std::vector<double> pz_vector(x + 2 * L, x + 3 * L);
//         std::vector<double> psi_vector(x + 3 * L, x + 4 * L);

//         std::vector<double> vx_vector(x + 4 * L, x + 5 * L);
//         std::vector<double> vy_vector(x + 5 * L, x + 6 * L);
//         std::vector<double> vz_vector(x + 6 * L, x + 7 * L);
//         std::vector<double> w_vector(x + 7 * L, x + 8 * L);

//         std::vector<double> ax_vector(x + 8  * L, x + 9  * L);
//         std::vector<double> ay_vector(x + 9  * L, x + 10 * L);
//         std::vector<double> az_vector(x + 10 * L, x + 11 * L);
//         std::vector<double> aw_vector(x + 11 * L, x + 12 * L);

//         std::vector<double> dyn1(L);
//         std::vector<double> dyn2(L);
//         std::vector<double> dyn3(L);
//         std::vector<double> dyn4(L);
//         std::vector<double> dyn5(L);
//         std::vector<double> dyn6(L);
//         std::vector<double> dyn7(L);
//         std::vector<double> dyn8(L);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, px_vector.data(), 1, 0.0, dyn1.data(), 1);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, py_vector.data(), 1, 0.0, dyn2.data(), 1);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, pz_vector.data(), 1, 0.0, dyn3.data(), 1);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, psi_vector.data(), 1, 0.0, dyn4.data(), 1);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, vx_vector.data(), 1, 0.0, dyn5.data(), 1);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, vy_vector.data(), 1, 0.0, dyn6.data(), 1);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, vz_vector.data(), 1, 0.0, dyn7.data(), 1);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, w_vector.data(), 1, 0.0, dyn8.data(), 1);

//         std::vector<double> px_rhs_scaled(L);
//         std::vector<double> py_rhs_scaled(L);
//         std::vector<double> pz_rhs_scaled(L);
//         std::vector<double> psi_rhs_scaled(L);

//         std::vector<double> vx_rhs_scaled(L);
//         std::vector<double> vy_rhs_scaled(L);
//         std::vector<double> vz_rhs_scaled(L);
//         std::vector<double> w_rhs_scaled(L);

//         vdMul(L, vx_vector.data(), r_tau_row_.data(), px_rhs_scaled.data());
//         vdMul(L, vy_vector.data(), r_tau_row_.data(), py_rhs_scaled.data());
//         vdMul(L, vz_vector.data(), r_tau_row_.data(), pz_rhs_scaled.data());
//         vdMul(L, w_vector.data(),  r_tau_row_.data(), psi_rhs_scaled.data());

//         vdMul(L, ax_vector.data(), r_tau_row_.data(), vx_rhs_scaled.data());
//         vdMul(L, ay_vector.data(), r_tau_row_.data(), vy_rhs_scaled.data());
//         vdMul(L, az_vector.data(), r_tau_row_.data(), vz_rhs_scaled.data());
//         vdMul(L, aw_vector.data(), r_tau_row_.data(), w_rhs_scaled.data());

//         std::vector<double> res_px(L);
//         std::vector<double> res_py(L);
//         std::vector<double> res_pz(L);
//         std::vector<double> res_psi(L);

//         std::vector<double> res_vx(L);
//         std::vector<double> res_vy(L);
//         std::vector<double> res_vz(L);
//         std::vector<double> res_w(L);

//         vdSub(L, dyn1.data(), px_rhs_scaled.data(),  res_px.data());
//         vdSub(L, dyn2.data(), py_rhs_scaled.data(),  res_py.data());
//         vdSub(L, dyn3.data(), pz_rhs_scaled.data(),  res_pz.data());
//         vdSub(L, dyn4.data(), psi_rhs_scaled.data(), res_psi.data());

//         vdSub(L, dyn5.data(), vx_rhs_scaled.data(), res_vx.data());
//         vdSub(L, dyn6.data(), vy_rhs_scaled.data(), res_vy.data());
//         vdSub(L, dyn7.data(), vz_rhs_scaled.data(), res_vz.data());
//         vdSub(L, dyn8.data(), w_rhs_scaled.data(),  res_w.data());

//         const double v_max2 = v_max_ * v_max_;
//         const double a_max2 = a_max_ * a_max_;

//         std::vector<double> vx2(L);
//         std::vector<double> vy2(L);
//         std::vector<double> vz2(L);

//         std::vector<double> ax2(L);
//         std::vector<double> ay2(L);
//         std::vector<double> az2(L);

//         vdSqr(L, vx_vector.data(), vx2.data());
//         vdSqr(L, vy_vector.data(), vy2.data());
//         vdSqr(L, vz_vector.data(), vz2.data());

//         vdSqr(L, ax_vector.data(), ax2.data());
//         vdSqr(L, ay_vector.data(), ay2.data());
//         vdSqr(L, az_vector.data(), az2.data());

//         std::vector<double> v2_sum(L);
//         std::vector<double> a2_sum(L);

//         vdAdd(L, vx2.data(), vy2.data(), v2_sum.data());
//         vdAdd(L, v2_sum.data(), vz2.data(), v2_sum.data());

//         vdAdd(L, ax2.data(), ay2.data(), a2_sum.data());
//         vdAdd(L, a2_sum.data(), az2.data(), a2_sum.data());

//         std::vector<double> v_max2_vec(L, v_max2);
//         std::vector<double> a_max2_vec(L, a_max2);

//         std::vector<double> c_speed(L);
//         std::vector<double> c_accel(L);

//         vdSub(L, v2_sum.data(), v_max2_vec.data(), c_speed.data());
//         vdSub(L, a2_sum.data(), a_max2_vec.data(), c_accel.data());

//         // ------------------------------------------------------------
//         // Obstacle avoidance constraints
//         //
//         // Same structure as MATLAB:
//         //
//         // dist2obs_square =
//         //     BernsteinProduct(px - ox, px - ox)
//         //   + BernsteinProduct(py - oy, py - oy)
//         //
//         // dist2obs_square_elev =
//         //     dist2obs_square * DegElevMatrix(...)
//         //
//         // c_obs = -dist2obs_square_elev + sep^2 <= 0
//         //
//         // Therefore:
//         //     dist2obs_square_elev >= sep^2
//         // ------------------------------------------------------------

//         auto add_vectors =
//             [](const std::vector<double>& a,
//             const std::vector<double>& b) -> std::vector<double> {
//                 const size_t size = std::min(a.size(), b.size());
//                 std::vector<double> out(size, 0.0);

//                 for (size_t i = 0; i < size; ++i) {
//                     out[i] = a[i] + b[i];
//                 }

//                 return out;
//             };

//         auto binom_ld =
//             [](int n, int k) -> long double {
//                 if (k < 0 || k > n) {
//                     return 0.0L;
//                 }

//                 if (k == 0 || k == n) {
//                     return 1.0L;
//                 }

//                 if (k > n - k) {
//                     k = n - k;
//                 }

//                 long double result = 1.0L;

//                 for (int i = 1; i <= k; ++i) {
//                     result *= static_cast<long double>(n - k + i);
//                     result /= static_cast<long double>(i);
//                 }

//                 return result;
//             };

//         auto degree_elevate =
//             [&](const std::vector<double>& cp,
//                 int target_degree) -> std::vector<double> {
//                 const int current_degree = static_cast<int>(cp.size()) - 1;

//                 if (target_degree <= current_degree) {
//                     return cp;
//                 }

//                 const int r = target_degree - current_degree;

//                 std::vector<double> elevated(target_degree + 1, 0.0);

//                 for (int j = 0; j <= target_degree; ++j) {
//                     const int i_min = std::max(0, j - r);
//                     const int i_max = std::min(current_degree, j);

//                     long double value = 0.0L;

//                     for (int i = i_min; i <= i_max; ++i) {
//                         const long double coeff =
//                             binom_ld(current_degree, i)
//                         * binom_ld(r, j - i)
//                         / binom_ld(target_degree, j);

//                         value += coeff * static_cast<long double>(cp[i]);
//                     }

//                     elevated[j] = static_cast<double>(value);
//                 }

//                 return elevated;
//             };

//         std::vector<double> c_obs_all;
//         const double sep2 = obs_sep_ * obs_sep_;

//         for (size_t obs_idx = 0; obs_idx < obs_x_.size(); ++obs_idx) {
//             const double ox = obs_x_[obs_idx];
//             const double oy = obs_y_[obs_idx];

//             std::vector<double> dx(L);
//             std::vector<double> dy(L);

//             for (int i = 0; i < L; ++i) {
//                 dx[i] = px_vector[i] - ox;
//                 dy[i] = py_vector[i] - oy;
//             }

//             std::vector<double> dx2 = BernsteinProduct(dx, dx);
//             std::vector<double> dy2 = BernsteinProduct(dy, dy);

//             std::vector<double> dist2obs_square =
//                 add_vectors(dx2, dy2);

//             const int deg_current =
//                 static_cast<int>(dist2obs_square.size()) - 1;

//             const int deg_target =
//                 deg_current + obs_deg_elev_extra_;

//             std::vector<double> dist2obs_square_elev =
//                 degree_elevate(dist2obs_square, deg_target);

//             for (double coeff : dist2obs_square_elev) {
//                 c_obs_all.push_back(-coeff + sep2);
//             }
//         }

//         for (Index i = 0; i < L; ++i) {
//             g[0 * L + i] = res_px[i];
//             g[1 * L + i] = res_py[i];
//             g[2 * L + i] = res_pz[i];
//             g[3 * L + i] = res_psi[i];

//             g[4 * L + i] = res_vx[i];
//             g[5 * L + i] = res_vy[i];
//             g[6 * L + i] = res_vz[i];
//             g[7 * L + i] = res_w[i];

//             g[8 * L + i] = c_speed[i];
//             g[9 * L + i] = c_accel[i];
//         }

//         const Index obs_start = 10 * L;

//         for (Index i = 0; i < static_cast<Index>(c_obs_all.size()); ++i) {
//             g[obs_start + i] = c_obs_all[i];
//         }

//         return true;
//     }

//     virtual bool eval_jac_g(
//         Index n,
//         const Number* x,
//         bool new_x,
//         Index m,
//         Index nele_jac,
//         Index* iRow,
//         Index* jCol,
//         Number* values
//     ) {
//         if (values == NULL) {
//             for (Index i = 0; i < m; i++) {
//                 for (Index j = 0; j < n; j++) {
//                     iRow[i * n + j] = i;
//                     jCol[i * n + j] = j;
//                 }
//             }
//         }

//         return true;
//     }

//     virtual bool eval_grad_f(
//         Index n,
//         const Number* x,
//         bool new_x,
//         Number* grad_f
//     ) {
//         return true;
//     }

//     virtual void finalize_solution(
//         SolverReturn status,
//         Index n,
//         const Number* x,
//         const Number* z_L,
//         const Number* z_U,
//         Index m,
//         const Number* g,
//         const Number* lambda,
//         Number obj_value,
//         const IpoptData* ip_data,
//         IpoptCalculatedQuantities* ip_cq
//     ) {
//         const int L = N_ + 1;

//         solution_x_.resize(12 * L);

//         for (Index i = 0; i < 12 * L; ++i) {
//             solution_x_[i] = x[i];
//         }

//         final_obj_value_ = obj_value;

//         /*
//             CHANGED ONLY THIS PART:

//             Before, finalize_solution() recalculated BeBOT:

//                 bebot_ = Bebot(N_, tf_);
//                 bebot_.calculate();

//             Since bebot_ was already calculated in the constructor, this
//             recalculation is removed.
//         */

//         const int K = 1000;

//         std::vector<double> t_norm(K);
//         std::vector<double> t_real(K);

//         for (int i = 0; i < K; ++i) {
//             const double s =
//                 static_cast<double>(i) / static_cast<double>(K - 1);

//             t_norm[i] = s;
//             t_real[i] = s * tf_;
//         }

//         std::vector<double> px_vector(
//             solution_x_.begin(),
//             solution_x_.begin() + L
//         );

//         std::vector<double> py_vector(
//             solution_x_.begin() + 1 * L,
//             solution_x_.begin() + 2 * L
//         );

//         std::vector<double> pz_vector(
//             solution_x_.begin() + 2 * L,
//             solution_x_.begin() + 3 * L
//         );

//         std::vector<double> psi_vector(
//             solution_x_.begin() + 3 * L,
//             solution_x_.begin() + 4 * L
//         );

//         std::vector<double> vx_vector(
//             solution_x_.begin() + 4 * L,
//             solution_x_.begin() + 5 * L
//         );

//         std::vector<double> vy_vector(
//             solution_x_.begin() + 5 * L,
//             solution_x_.begin() + 6 * L
//         );

//         std::vector<double> vz_vector(
//             solution_x_.begin() + 6 * L,
//             solution_x_.begin() + 7 * L
//         );

//         std::vector<double> w_vector(
//             solution_x_.begin() + 7 * L,
//             solution_x_.begin() + 8 * L
//         );

//         std::vector<double> ax_vector(
//             solution_x_.begin() + 8 * L,
//             solution_x_.begin() + 9 * L
//         );

//         std::vector<double> ay_vector(
//             solution_x_.begin() + 9 * L,
//             solution_x_.begin() + 10 * L
//         );

//         std::vector<double> az_vector(
//             solution_x_.begin() + 10 * L,
//             solution_x_.begin() + 11 * L
//         );

//         std::vector<double> aw_vector(
//             solution_x_.begin() + 11 * L,
//             solution_x_.begin() + 12 * L
//         );

//         std::vector<std::vector<double>> px_2d(1, px_vector);
//         std::vector<std::vector<double>> py_2d(1, py_vector);
//         std::vector<std::vector<double>> pz_2d(1, pz_vector);
//         std::vector<std::vector<double>> psi_2d(1, psi_vector);

//         std::vector<std::vector<double>> vx_2d(1, vx_vector);
//         std::vector<std::vector<double>> vy_2d(1, vy_vector);
//         std::vector<std::vector<double>> vz_2d(1, vz_vector);
//         std::vector<std::vector<double>> w_2d(1, w_vector);

//         std::vector<std::vector<double>> ax_2d(1, ax_vector);
//         std::vector<std::vector<double>> ay_2d(1, ay_vector);
//         std::vector<std::vector<double>> az_2d(1, az_vector);
//         std::vector<std::vector<double>> aw_2d(1, aw_vector);

//         auto px_real =
//             BernsteinPoly(px_2d, t_real, 0.0, tf_);

//         auto py_real =
//             BernsteinPoly(py_2d, t_real, 0.0, tf_);

//         auto pz_real =
//             BernsteinPoly(pz_2d, t_real, 0.0, tf_);

//         auto psi_real =
//             BernsteinPoly(psi_2d, t_real, 0.0, tf_);

//         auto vx_real =
//             BernsteinPoly(vx_2d, t_real, 0.0, tf_);

//         auto vy_real =
//             BernsteinPoly(vy_2d, t_real, 0.0, tf_);

//         auto vz_real =
//             BernsteinPoly(vz_2d, t_real, 0.0, tf_);

//         auto w_real =
//             BernsteinPoly(w_2d, t_real, 0.0, tf_);

//         auto ax_real =
//             BernsteinPoly(ax_2d, t_real, 0.0, tf_);

//         auto ay_real =
//             BernsteinPoly(ay_2d, t_real, 0.0, tf_);

//         auto az_real =
//             BernsteinPoly(az_2d, t_real, 0.0, tf_);

//         auto aw_real =
//             BernsteinPoly(aw_2d, t_real, 0.0, tf_);

//         auto flatten =
//             [](const std::vector<std::vector<double>>& input) {
//                 std::vector<double> output;

//                 for (const auto& row : input) {
//                     output.insert(output.end(), row.begin(), row.end());
//                 }

//                 return output;
//             };

//         writeToCSV(t_real, flatten(px_real), "px_real.csv");
//         writeToCSV(bebot_.getNodes(), px_vector, "px_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(py_real), "py_real.csv");
//         writeToCSV(bebot_.getNodes(), py_vector, "py_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(pz_real), "pz_real.csv");
//         writeToCSV(bebot_.getNodes(), pz_vector, "pz_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(psi_real), "psi_real.csv");
//         writeToCSV(bebot_.getNodes(), psi_vector, "psi_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(vx_real), "vx_real.csv");
//         writeToCSV(bebot_.getNodes(), vx_vector, "vx_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(vy_real), "vy_real.csv");
//         writeToCSV(bebot_.getNodes(), vy_vector, "vy_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(vz_real), "vz_real.csv");
//         writeToCSV(bebot_.getNodes(), vz_vector, "vz_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(w_real), "w_real.csv");
//         writeToCSV(bebot_.getNodes(), w_vector, "w_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(ax_real), "ax_real.csv");
//         writeToCSV(bebot_.getNodes(), ax_vector, "ax_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(ay_real), "ay_real.csv");
//         writeToCSV(bebot_.getNodes(), ay_vector, "ay_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(az_real), "az_real.csv");
//         writeToCSV(bebot_.getNodes(), az_vector, "az_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(aw_real), "aw_real.csv");
//         writeToCSV(bebot_.getNodes(), aw_vector, "aw_controlpoints_real.csv");
//     }

//     const std::vector<Number>& get_solution_x() const {
//         return solution_x_;
//     }

//     Number get_final_obj_value() const {
//         return final_obj_value_;
//     }

// private:
//     int N_;
//     double tf_;

//     double px_max_;
//     double px_min_;
//     double py_max_;
//     double py_min_;
//     double pz_max_;
//     double pz_min_;
//     double psi_max_;
//     double psi_min_;

//     double v_max_;
//     double w_max_;
//     double a_max_;
//     double aw_max_;

//     double px_cur_;
//     double py_cur_;
//     double pz_cur_;
//     double psi_cur_;

//     double vx_cur_;
//     double vy_cur_;
//     double vz_cur_;
//     double w_cur_;

//     double pxf_;
//     double pyf_;
//     double pzf_;
//     double psif_;

//     std::array<double, 2> obs_x_;
//     std::array<double, 2> obs_y_;
//     double obs_sep_;
//     int obs_deg_elev_extra_;

//     Bebot bebot_;

//     std::vector<double> tau_vec_;
//     std::vector<double> r_tau_;
//     std::vector<double> r_tau_row_;

//     std::vector<Number> solution_u_;
//     std::vector<Number> solution_x2_;
//     std::vector<Number> solution_x_;

//     Number final_obj_value_;

//     std::vector<double> final_time_;

//     std::vector<std::vector<double>> bernsteinpoly_resultu_;
//     std::vector<std::vector<double>> bernsteinpoly_resultx2_;
//     std::vector<std::vector<double>> bernsteinpoly_resultz_;

// public:
//     const std::vector<std::vector<double>>& get_bernsteinpoly_result() const {
//         return bernsteinpoly_resultz_;
//     }
// };

// extern "C" {
//     PointSetProblem* create_point_set_problem(
//         int N,
//         double tf,
//         double px_max,
//         double px_min,
//         double py_max,
//         double py_min,
//         double pz_max,
//         double pz_min,
//         double psi_max,
//         double psi_min,
//         double v_max,
//         double w_max,
//         double a_max,
//         double aw_max,
//         double px_cur,
//         double py_cur,
//         double pz_cur,
//         double psi_cur,
//         double vx_cur,
//         double vy_cur,
//         double vz_cur,
//         double w_cur,
//         double pxf,
//         double pyf,
//         double pzf,
//         double psif
//     ) {
//         return new PointSetProblem(
//             N,
//             tf,
//             px_max,
//             px_min,
//             py_max,
//             py_min,
//             pz_max,
//             pz_min,
//             psi_max,
//             psi_min,
//             v_max,
//             w_max,
//             a_max,
//             aw_max,
//             px_cur,
//             py_cur,
//             pz_cur,
//             psi_cur,
//             vx_cur,
//             vy_cur,
//             vz_cur,
//             w_cur,
//             pxf,
//             pyf,
//             pzf,
//             psif
//         );
//     }

//     void solve_point_set_problem(PointSetProblem* problem) {
//         SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

//         app->Options()->SetStringValue("linear_solver", "ma57");
//         app->Options()->SetStringValue("mu_strategy", "adaptive");

//         app->Options()->SetStringValue(
//             "gradient_approximation",
//             "finite-difference-values"
//         );

//         app->Options()->SetStringValue(
//             "jacobian_approximation",
//             "finite-difference-values"
//         );

//         app->Options()->SetStringValue(
//             "hessian_approximation",
//             "limited-memory"
//         );

//         app->Options()->SetIntegerValue("max_iter", 400);

//         app->Options()->SetNumericValue("tol", 1e-3);
//         app->Options()->SetNumericValue("constr_viol_tol", 1e-3);
//         app->Options()->SetNumericValue("obj_scaling_factor", 1e-3);

//         app->Options()->SetIntegerValue("print_level", 3);

//         app->RethrowNonIpoptException(true);

//         ApplicationReturnStatus status = app->Initialize();

//         if (status != Solve_Succeeded) {
//             std::cerr << "IPOPT initialization failed!" << std::endl;
//             return;
//         }

//         status = app->OptimizeTNLP(problem);

//         if (status == Solve_Succeeded || status == Solved_To_Acceptable_Level) {
//             std::cout << "Optimization succeeded!" << std::endl;

//             const auto& solution_x = problem->get_solution_x();

//             std::cout << "Optimal Solution (x): ";

//             for (Index i = 0; i < static_cast<Index>(solution_x.size()); i++) {
//                 std::cout << solution_x[i] << " ";
//             }

//             std::cout << std::endl;
//         } else {
//             std::cerr << "Optimization failed with status "
//                       << status << std::endl;
//         }
//     }

//     void get_solution(PointSetProblem* problem, double* solution, int n) {
//         const std::vector<double>& sol = problem->get_solution_x();

//         std::cout << "[DEBUG] get_solution called. Vector size: "
//                   << sol.size() << std::endl;

//         std::copy(sol.begin(), sol.end(), solution);
//     }

//     double get_final_objective_value(PointSetProblem* problem) {
//         return problem->get_final_obj_value();
//     }

//     void destroy_point_set_problem(PointSetProblem* problem) {
//         delete problem;
//     }
// }

////////////////////////////////////////////////////////////////////////////////////////////////////////

// #include "../../../../Ipopt_ma57_solver/src/Interfaces/IpIpoptApplication.hpp"
// #include "../../../../Ipopt_ma57_solver/src/Interfaces/IpTNLP.hpp"

// #include <cmath>
// #include <iostream>
// #include <fstream>
// #include <vector>
// #include <iomanip>
// #include <array>
// #include <string>
// #include <limits>
// #include <algorithm>
// #include <sstream>

// #include "../../../../include/bebot.h"
// #include "../../../../include/bernsteinpoly.h"
// #include "../../../../include/bernsteinproduct.h"
// #include "../../../../include/degelevmatrix.h"

// #include "mkl.h"

// using namespace Ipopt;

// class PointSetProblem : public Ipopt::TNLP {
// public:
//     PointSetProblem(
//         int N,
//         double tf,
//         double px_max,
//         double px_min,
//         double py_max,
//         double py_min,
//         double pz_max,
//         double pz_min,
//         double psi_max,
//         double psi_min,
//         double v_max,
//         double w_max,
//         double a_max,
//         double aw_max,
//         double px_cur,
//         double py_cur,
//         double pz_cur,
//         double psi_cur,
//         double vx_cur,
//         double vy_cur,
//         double vz_cur,
//         double w_cur,
//         double pxf,
//         double pyf,
//         double pzf,
//         double psif
//     )
//         : N_(N),
//           tf_(tf),
//           px_max_(px_max),
//           px_min_(px_min),
//           py_max_(py_max),
//           py_min_(py_min),
//           pz_max_(pz_max),
//           pz_min_(pz_min),
//           psi_max_(psi_max),
//           psi_min_(psi_min),
//           v_max_(v_max),
//           w_max_(w_max),
//           a_max_(a_max),
//           aw_max_(aw_max),
//           px_cur_(px_cur),
//           py_cur_(py_cur),
//           pz_cur_(pz_cur),
//           psi_cur_(psi_cur),
//           vx_cur_(vx_cur),
//           vy_cur_(vy_cur),
//           vz_cur_(vz_cur),
//           w_cur_(w_cur),
//           pxf_(pxf),
//           pyf_(pyf),
//           pzf_(pzf),
//           psif_(psif),
//           obs_x_{{0.750}},
//           obs_y_{{ 1.4}},
//           obs_sep_(0.50),
//           obs_deg_elev_extra_(5),
//           bebot_(N, tf_) {
//         bebot_.calculate();

//         const int L = N_ + 1;

//         tau_vec_.resize(L);
//         r_tau_.resize(L);
//         r_tau_row_.resize(L);

//         if (L == 1) {
//             const double tau0 = -1.0;
//             const double denom = 1.0 - tau0;

//             tau_vec_[0] = tau0;
//             r_tau_[0] = 2.0 / (denom * denom);
//             r_tau_row_[0] = r_tau_[0];
//         } else {
//             const double tau0 = -1.0;
//             const double tau1 = 0.9;
//             const double step = (tau1 - tau0) / static_cast<double>(N_);

//             for (int i = 0; i < L; ++i) {
//                 const double tau = tau0 + step * static_cast<double>(i);
//                 const double denom = 1.0 - tau;

//                 tau_vec_[i] = tau;
//                 r_tau_[i] = 2.0 / (denom * denom);
//                 r_tau_row_[i] = r_tau_[i];
//             }
//         }
//     }

//     void writeToCSV(
//         const std::vector<double>& times,
//         const std::vector<double>& values,
//         const std::string& filename
//     ) {
//         std::ofstream outFile(filename);

//         if (!outFile.is_open()) {
//             std::cerr << "Failed to open file: " << filename << std::endl;
//             return;
//         }

//         outFile << "Time,Value\n";

//         for (size_t i = 0; i < times.size(); ++i) {
//             outFile << std::fixed << std::setprecision(6)
//                     << times[i] << "," << values[i] << "\n";
//         }

//         outFile.close();
//     }

//     virtual bool get_nlp_info(
//         Index& n,
//         Index& m,
//         Index& nnz_jac_g,
//         Index& nnz_h_lag,
//         IndexStyleEnum& index_style
//     ) {
//         const int L = N_ + 1;

//         n = 12 * L;

//         /*
//             Existing constraints:
//             8 * (N + 1): dynamics equality constraints
//             1 * (N + 1): speed inequality constraints
//             1 * (N + 1): acceleration inequality constraints

//             Obstacle constraints:
//             For each obstacle:
//                 dist2obs_square has degree 2N
//                 then degree-elevated to 2N + obs_deg_elev_extra_
//                 number of coefficients = 2N + obs_deg_elev_extra_ + 1

//             Total:
//                 m = 10 * (N + 1)
//                     + number_of_obstacles * (2N + obs_deg_elev_extra_ + 1)
//         */

//         const int obs_degree = 2 * N_ + obs_deg_elev_extra_;
//         const int obs_constraints =
//             static_cast<int>(obs_x_.size()) * (obs_degree + 1);

//         m = 10 * L + obs_constraints;

//         nnz_jac_g = n * m;
//         nnz_h_lag = 0;
//         index_style = TNLP::C_STYLE;

//         return true;
//     }

//     virtual bool get_bounds_info(
//         Index n,
//         Number* x_l,
//         Number* x_u,
//         Index m,
//         Number* g_l,
//         Number* g_u
//     ) {
//         const int L = N_ + 1;

//         std::vector<double> x_lower(
//             n,
//             -std::numeric_limits<double>::infinity()
//         );

//         std::vector<double> x_upper(
//             n,
//             std::numeric_limits<double>::infinity()
//         );

//         // px
//         for (int i = 1; i < L; ++i) {
//             x_lower[i] = px_min_;
//             x_upper[i] = px_max_;
//         }
//         x_lower[0] = px_cur_;
//         x_upper[0] = px_cur_;

//         // py
//         for (int i = L + 1; i < 2 * L; ++i) {
//             x_lower[i] = py_min_;
//             x_upper[i] = py_max_;
//         }
//         x_lower[L] = py_cur_;
//         x_upper[L] = py_cur_;

//         // pz
//         for (int i = 2 * L + 1; i < 3 * L; ++i) {
//             x_lower[i] = pz_min_;
//             x_upper[i] = pz_max_;
//         }
//         x_lower[2 * L] = pz_cur_;
//         x_upper[2 * L] = pz_cur_;

//         // psi
//         for (int i = 3 * L + 1; i < 4 * L; ++i) {
//             x_lower[i] = psi_min_;
//             x_upper[i] = psi_max_;
//         }
//         x_lower[3 * L] = psi_cur_;
//         x_upper[3 * L] = psi_cur_;

//         // vx
//         for (int i = 4 * L + 1; i < 5 * L; ++i) {
//             x_lower[i] = -v_max_;
//             x_upper[i] =  v_max_;
//         }
//         x_lower[4 * L] = vx_cur_;
//         x_upper[4 * L] = vx_cur_;

//         // vy
//         for (int i = 5 * L + 1; i < 6 * L; ++i) {
//             x_lower[i] = -v_max_;
//             x_upper[i] =  v_max_;
//         }
//         x_lower[5 * L] = vy_cur_;
//         x_upper[5 * L] = vy_cur_;

//         // vz
//         for (int i = 6 * L + 1; i < 7 * L; ++i) {
//             x_lower[i] = -v_max_;
//             x_upper[i] =  v_max_;
//         }
//         x_lower[6 * L] = vz_cur_;
//         x_upper[6 * L] = vz_cur_;

//         // w
//         for (int i = 7 * L + 1; i < 8 * L; ++i) {
//             x_lower[i] = -w_max_;
//             x_upper[i] =  w_max_;
//         }
//         x_lower[7 * L] = w_cur_;
//         x_upper[7 * L] = w_cur_;

//         // ax
//         for (int i = 8 * L; i < 9 * L; ++i) {
//             x_lower[i] = -a_max_;
//             x_upper[i] =  a_max_;
//         }

//         // ay
//         for (int i = 9 * L; i < 10 * L; ++i) {
//             x_lower[i] = -a_max_;
//             x_upper[i] =  a_max_;
//         }

//         // az
//         for (int i = 10 * L; i < 11 * L; ++i) {
//             x_lower[i] = -a_max_;
//             x_upper[i] =  a_max_;
//         }

//         // aw
//         for (int i = 11 * L; i < 12 * L; ++i) {
//             x_lower[i] = -aw_max_;
//             x_upper[i] =  aw_max_;
//         }

//         std::copy(x_lower.begin(), x_lower.end(), x_l);
//         std::copy(x_upper.begin(), x_upper.end(), x_u);

//         // Dynamics equality constraints
//         for (int i = 0; i < 8 * L; ++i) {
//             g_l[i] = 0.0;
//             g_u[i] = 0.0;
//         }

//         // Speed inequality constraints: c_speed <= 0
//         for (int i = 8 * L; i < 9 * L; ++i) {
//             g_l[i] = -std::numeric_limits<double>::infinity();
//             g_u[i] = 0.0;
//         }

//         // Acceleration inequality constraints: c_accel <= 0
//         for (int i = 9 * L; i < 10 * L; ++i) {
//             g_l[i] = -std::numeric_limits<double>::infinity();
//             g_u[i] = 0.0;
//         }

//         // Obstacle inequality constraints: c_obs <= 0
//         const int obs_start = 10 * L;

//         for (Index i = obs_start; i < m; ++i) {
//             g_l[i] = -std::numeric_limits<double>::infinity();
//             g_u[i] = 0.0;
//         }

//         return true;
//     }

//     virtual bool get_starting_point(
//         Index n,
//         bool init_x,
//         Number* x,
//         bool init_z,
//         Number* z_L,
//         Number* z_U,
//         Index m,
//         bool init_lambda,
//         Number* lambda
//     ) {
//         for (Index i = 0; i < n; ++i) {
//             x[i] = 1.0;
//         }

//         return true;
//     }

//     virtual bool eval_f(
//         Index n,
//         const Number* x,
//         bool new_x,
//         Number& obj_value
//     ) {
//         const int L = N_ + 1;

//         const double w_p   = 0.01;
//         const double w_psi = 1.0;
//         const double w_a   = 0.01;
//         const double w_aw  = 0.001;

//         obj_value = 0.0;

//         std::vector<double> pxf_vector(L, pxf_);
//         std::vector<double> pyf_vector(L, pyf_);
//         std::vector<double> pzf_vector(L, pzf_);
//         std::vector<double> psif_vector(L, psif_);

//         std::vector<double> px_vector(x, x + L);
//         std::vector<double> py_vector(x + 1 * L, x + 2 * L);
//         std::vector<double> pz_vector(x + 2 * L, x + 3 * L);
//         std::vector<double> psi_vector(x + 3 * L, x + 4 * L);

//         std::vector<double> ax_vector(x + 8  * L, x + 9  * L);
//         std::vector<double> ay_vector(x + 9  * L, x + 10 * L);
//         std::vector<double> az_vector(x + 10 * L, x + 11 * L);
//         std::vector<double> aw_vector(x + 11 * L, x + 12 * L);

//         std::vector<double> px_diff(L);
//         std::vector<double> py_diff(L);
//         std::vector<double> pz_diff(L);
//         std::vector<double> psi_diff(L);

//         vdSub(L, px_vector.data(),  pxf_vector.data(),  px_diff.data());
//         vdSub(L, py_vector.data(),  pyf_vector.data(),  py_diff.data());
//         vdSub(L, pz_vector.data(),  pzf_vector.data(),  pz_diff.data());
//         vdSub(L, psi_vector.data(), psif_vector.data(), psi_diff.data());

//         std::vector<double> px_diff_sqr(L);
//         std::vector<double> py_diff_sqr(L);
//         std::vector<double> pz_diff_sqr(L);
//         std::vector<double> psi_diff_sqr(L);

//         vdSqr(L, px_diff.data(),  px_diff_sqr.data());
//         vdSqr(L, py_diff.data(),  py_diff_sqr.data());
//         vdSqr(L, pz_diff.data(),  pz_diff_sqr.data());
//         vdSqr(L, psi_diff.data(), psi_diff_sqr.data());

//         std::vector<double> px_weighted(L);
//         std::vector<double> py_weighted(L);
//         std::vector<double> pz_weighted(L);
//         std::vector<double> psi_weighted(L);

//         vdMul(L, r_tau_.data(), px_diff_sqr.data(),  px_weighted.data());
//         vdMul(L, r_tau_.data(), py_diff_sqr.data(),  py_weighted.data());
//         vdMul(L, r_tau_.data(), pz_diff_sqr.data(),  pz_weighted.data());
//         vdMul(L, r_tau_.data(), psi_diff_sqr.data(), psi_weighted.data());

//         std::vector<double> ones(L, 1.0);

//         const double sum_px_state =
//             cblas_ddot(L, px_weighted.data(), 1, ones.data(), 1);

//         const double sum_py_state =
//             cblas_ddot(L, py_weighted.data(), 1, ones.data(), 1);

//         const double sum_pz_state =
//             cblas_ddot(L, pz_weighted.data(), 1, ones.data(), 1);

//         const double sum_psi_state =
//             cblas_ddot(L, psi_weighted.data(), 1, ones.data(), 1);

//         const double state_term =
//             w_p   * sum_px_state
//           + w_p   * sum_py_state
//           + w_p   * sum_pz_state
//           + w_psi * sum_psi_state;

//         std::vector<double> ax_sqr(L);
//         std::vector<double> ay_sqr(L);
//         std::vector<double> az_sqr(L);
//         std::vector<double> aw_sqr(L);

//         vdSqr(L, ax_vector.data(), ax_sqr.data());
//         vdSqr(L, ay_vector.data(), ay_sqr.data());
//         vdSqr(L, az_vector.data(), az_sqr.data());
//         vdSqr(L, aw_vector.data(), aw_sqr.data());

//         const double sum_ax =
//             cblas_ddot(L, ax_sqr.data(), 1, ones.data(), 1);

//         const double sum_ay =
//             cblas_ddot(L, ay_sqr.data(), 1, ones.data(), 1);

//         const double sum_az =
//             cblas_ddot(L, az_sqr.data(), 1, ones.data(), 1);

//         const double sum_aw =
//             cblas_ddot(L, aw_sqr.data(), 1, ones.data(), 1);

//         const double control_term =
//             w_a  * sum_ax
//           + w_a  * sum_ay
//           + w_a  * sum_az
//           + w_aw * sum_aw;

//         obj_value = state_term + control_term;

//         return true;
//     }

//     virtual bool eval_g(
//         Index n,
//         const Number* x,
//         bool new_x,
//         Index m,
//         Number* g
//     ) {
//         const int L = N_ + 1;

//         const auto& Dm = bebot_.getDifferentiationMatrix();

//         std::vector<double> px_vector(x, x + L);
//         std::vector<double> py_vector(x + 1 * L, x + 2 * L);
//         std::vector<double> pz_vector(x + 2 * L, x + 3 * L);
//         std::vector<double> psi_vector(x + 3 * L, x + 4 * L);

//         std::vector<double> vx_vector(x + 4 * L, x + 5 * L);
//         std::vector<double> vy_vector(x + 5 * L, x + 6 * L);
//         std::vector<double> vz_vector(x + 6 * L, x + 7 * L);
//         std::vector<double> w_vector(x + 7 * L, x + 8 * L);

//         std::vector<double> ax_vector(x + 8  * L, x + 9  * L);
//         std::vector<double> ay_vector(x + 9  * L, x + 10 * L);
//         std::vector<double> az_vector(x + 10 * L, x + 11 * L);
//         std::vector<double> aw_vector(x + 11 * L, x + 12 * L);

//         std::vector<double> dyn1(L);
//         std::vector<double> dyn2(L);
//         std::vector<double> dyn3(L);
//         std::vector<double> dyn4(L);
//         std::vector<double> dyn5(L);
//         std::vector<double> dyn6(L);
//         std::vector<double> dyn7(L);
//         std::vector<double> dyn8(L);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, px_vector.data(), 1, 0.0, dyn1.data(), 1);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, py_vector.data(), 1, 0.0, dyn2.data(), 1);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, pz_vector.data(), 1, 0.0, dyn3.data(), 1);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, psi_vector.data(), 1, 0.0, dyn4.data(), 1);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, vx_vector.data(), 1, 0.0, dyn5.data(), 1);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, vy_vector.data(), 1, 0.0, dyn6.data(), 1);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, vz_vector.data(), 1, 0.0, dyn7.data(), 1);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, w_vector.data(), 1, 0.0, dyn8.data(), 1);

//         std::vector<double> px_rhs_scaled(L);
//         std::vector<double> py_rhs_scaled(L);
//         std::vector<double> pz_rhs_scaled(L);
//         std::vector<double> psi_rhs_scaled(L);

//         std::vector<double> vx_rhs_scaled(L);
//         std::vector<double> vy_rhs_scaled(L);
//         std::vector<double> vz_rhs_scaled(L);
//         std::vector<double> w_rhs_scaled(L);

//         vdMul(L, vx_vector.data(), r_tau_row_.data(), px_rhs_scaled.data());
//         vdMul(L, vy_vector.data(), r_tau_row_.data(), py_rhs_scaled.data());
//         vdMul(L, vz_vector.data(), r_tau_row_.data(), pz_rhs_scaled.data());
//         vdMul(L, w_vector.data(),  r_tau_row_.data(), psi_rhs_scaled.data());

//         vdMul(L, ax_vector.data(), r_tau_row_.data(), vx_rhs_scaled.data());
//         vdMul(L, ay_vector.data(), r_tau_row_.data(), vy_rhs_scaled.data());
//         vdMul(L, az_vector.data(), r_tau_row_.data(), vz_rhs_scaled.data());
//         vdMul(L, aw_vector.data(), r_tau_row_.data(), w_rhs_scaled.data());

//         std::vector<double> res_px(L);
//         std::vector<double> res_py(L);
//         std::vector<double> res_pz(L);
//         std::vector<double> res_psi(L);

//         std::vector<double> res_vx(L);
//         std::vector<double> res_vy(L);
//         std::vector<double> res_vz(L);
//         std::vector<double> res_w(L);

//         vdSub(L, dyn1.data(), px_rhs_scaled.data(),  res_px.data());
//         vdSub(L, dyn2.data(), py_rhs_scaled.data(),  res_py.data());
//         vdSub(L, dyn3.data(), pz_rhs_scaled.data(),  res_pz.data());
//         vdSub(L, dyn4.data(), psi_rhs_scaled.data(), res_psi.data());

//         vdSub(L, dyn5.data(), vx_rhs_scaled.data(), res_vx.data());
//         vdSub(L, dyn6.data(), vy_rhs_scaled.data(), res_vy.data());
//         vdSub(L, dyn7.data(), vz_rhs_scaled.data(), res_vz.data());
//         vdSub(L, dyn8.data(), w_rhs_scaled.data(),  res_w.data());

//         const double v_max2 = v_max_ * v_max_;
//         const double a_max2 = a_max_ * a_max_;

//         std::vector<double> vx2(L);
//         std::vector<double> vy2(L);
//         std::vector<double> vz2(L);

//         std::vector<double> ax2(L);
//         std::vector<double> ay2(L);
//         std::vector<double> az2(L);

//         vdSqr(L, vx_vector.data(), vx2.data());
//         vdSqr(L, vy_vector.data(), vy2.data());
//         vdSqr(L, vz_vector.data(), vz2.data());

//         vdSqr(L, ax_vector.data(), ax2.data());
//         vdSqr(L, ay_vector.data(), ay2.data());
//         vdSqr(L, az_vector.data(), az2.data());

//         std::vector<double> v2_sum(L);
//         std::vector<double> a2_sum(L);

//         vdAdd(L, vx2.data(), vy2.data(), v2_sum.data());
//         vdAdd(L, v2_sum.data(), vz2.data(), v2_sum.data());

//         vdAdd(L, ax2.data(), ay2.data(), a2_sum.data());
//         vdAdd(L, a2_sum.data(), az2.data(), a2_sum.data());

//         std::vector<double> v_max2_vec(L, v_max2);
//         std::vector<double> a_max2_vec(L, a_max2);

//         std::vector<double> c_speed(L);
//         std::vector<double> c_accel(L);

//         vdSub(L, v2_sum.data(), v_max2_vec.data(), c_speed.data());
//         vdSub(L, a2_sum.data(), a_max2_vec.data(), c_accel.data());

//         auto add_vectors =
//             [](const std::vector<double>& a,
//             const std::vector<double>& b) -> std::vector<double> {
//                 const size_t size = std::min(a.size(), b.size());
//                 std::vector<double> out(size, 0.0);

//                 for (size_t i = 0; i < size; ++i) {
//                     out[i] = a[i] + b[i];
//                 }

//                 return out;
//             };

//         auto binom_ld =
//             [](int n, int k) -> long double {
//                 if (k < 0 || k > n) {
//                     return 0.0L;
//                 }

//                 if (k == 0 || k == n) {
//                     return 1.0L;
//                 }

//                 if (k > n - k) {
//                     k = n - k;
//                 }

//                 long double result = 1.0L;

//                 for (int i = 1; i <= k; ++i) {
//                     result *= static_cast<long double>(n - k + i);
//                     result /= static_cast<long double>(i);
//                 }

//                 return result;
//             };

//         auto degree_elevate =
//             [&](const std::vector<double>& cp,
//                 int target_degree) -> std::vector<double> {
//                 const int current_degree = static_cast<int>(cp.size()) - 1;

//                 if (target_degree <= current_degree) {
//                     return cp;
//                 }

//                 const int r = target_degree - current_degree;

//                 std::vector<double> elevated(target_degree + 1, 0.0);

//                 for (int j = 0; j <= target_degree; ++j) {
//                     const int i_min = std::max(0, j - r);
//                     const int i_max = std::min(current_degree, j);

//                     long double value = 0.0L;

//                     for (int i = i_min; i <= i_max; ++i) {
//                         const long double coeff =
//                             binom_ld(current_degree, i)
//                         * binom_ld(r, j - i)
//                         / binom_ld(target_degree, j);

//                         value += coeff * static_cast<long double>(cp[i]);
//                     }

//                     elevated[j] = static_cast<double>(value);
//                 }

//                 return elevated;
//             };

//         std::vector<double> c_obs_all;
//         const double sep2 = obs_sep_ * obs_sep_;

//         for (size_t obs_idx = 0; obs_idx < obs_x_.size(); ++obs_idx) {
//             const double ox = obs_x_[obs_idx];
//             const double oy = obs_y_[obs_idx];

//             std::vector<double> dx(L);
//             std::vector<double> dy(L);

//             for (int i = 0; i < L; ++i) {
//                 dx[i] = px_vector[i] - ox;
//                 dy[i] = py_vector[i] - oy;
//             }

//             std::vector<double> dx2 = BernsteinProduct(dx, dx);
//             std::vector<double> dy2 = BernsteinProduct(dy, dy);

//             std::vector<double> dist2obs_square =
//                 add_vectors(dx2, dy2);

//             const int deg_current =
//                 static_cast<int>(dist2obs_square.size()) - 1;

//             const int deg_target =
//                 deg_current + obs_deg_elev_extra_;

//             std::vector<double> dist2obs_square_elev =
//                 degree_elevate(dist2obs_square, deg_target);

//             for (double coeff : dist2obs_square_elev) {
//                 c_obs_all.push_back(-coeff + sep2);
//             }
//         }

//         for (Index i = 0; i < L; ++i) {
//             g[0 * L + i] = res_px[i];
//             g[1 * L + i] = res_py[i];
//             g[2 * L + i] = res_pz[i];
//             g[3 * L + i] = res_psi[i];

//             g[4 * L + i] = res_vx[i];
//             g[5 * L + i] = res_vy[i];
//             g[6 * L + i] = res_vz[i];
//             g[7 * L + i] = res_w[i];

//             g[8 * L + i] = c_speed[i];
//             g[9 * L + i] = c_accel[i];
//         }

//         const Index obs_start = 10 * L;

//         for (Index i = 0; i < static_cast<Index>(c_obs_all.size()); ++i) {
//             g[obs_start + i] = c_obs_all[i];
//         }

//         return true;
//     }

//     virtual bool eval_jac_g(
//         Index n,
//         const Number* x,
//         bool new_x,
//         Index m,
//         Index nele_jac,
//         Index* iRow,
//         Index* jCol,
//         Number* values
//     ) {
//         if (values == NULL) {
//             for (Index i = 0; i < m; i++) {
//                 for (Index j = 0; j < n; j++) {
//                     iRow[i * n + j] = i;
//                     jCol[i * n + j] = j;
//                 }
//             }
//         }

//         return true;
//     }

//     virtual bool eval_grad_f(
//         Index n,
//         const Number* x,
//         bool new_x,
//         Number* grad_f
//     ) {
//         return true;
//     }

//     virtual void finalize_solution(
//         SolverReturn status,
//         Index n,
//         const Number* x,
//         const Number* z_L,
//         const Number* z_U,
//         Index m,
//         const Number* g,
//         const Number* lambda,
//         Number obj_value,
//         const IpoptData* ip_data,
//         IpoptCalculatedQuantities* ip_cq
//     ) {
//         const int L = N_ + 1;

//         solution_x_.resize(12 * L);

//         for (Index i = 0; i < 12 * L; ++i) {
//             solution_x_[i] = x[i];
//         }

//         final_obj_value_ = obj_value;

//         const int K = 1000;

//         std::vector<double> t_norm(K);
//         std::vector<double> t_real(K);

//         for (int i = 0; i < K; ++i) {
//             const double s =
//                 static_cast<double>(i) / static_cast<double>(K - 1);

//             t_norm[i] = s;
//             t_real[i] = s * tf_;
//         }

//         std::vector<double> px_vector(
//             solution_x_.begin(),
//             solution_x_.begin() + L
//         );

//         std::vector<double> py_vector(
//             solution_x_.begin() + 1 * L,
//             solution_x_.begin() + 2 * L
//         );

//         std::vector<double> pz_vector(
//             solution_x_.begin() + 2 * L,
//             solution_x_.begin() + 3 * L
//         );

//         std::vector<double> psi_vector(
//             solution_x_.begin() + 3 * L,
//             solution_x_.begin() + 4 * L
//         );

//         std::vector<double> vx_vector(
//             solution_x_.begin() + 4 * L,
//             solution_x_.begin() + 5 * L
//         );

//         std::vector<double> vy_vector(
//             solution_x_.begin() + 5 * L,
//             solution_x_.begin() + 6 * L
//         );

//         std::vector<double> vz_vector(
//             solution_x_.begin() + 6 * L,
//             solution_x_.begin() + 7 * L
//         );

//         std::vector<double> w_vector(
//             solution_x_.begin() + 7 * L,
//             solution_x_.begin() + 8 * L
//         );

//         std::vector<double> ax_vector(
//             solution_x_.begin() + 8 * L,
//             solution_x_.begin() + 9 * L
//         );

//         std::vector<double> ay_vector(
//             solution_x_.begin() + 9 * L,
//             solution_x_.begin() + 10 * L
//         );

//         std::vector<double> az_vector(
//             solution_x_.begin() + 10 * L,
//             solution_x_.begin() + 11 * L
//         );

//         std::vector<double> aw_vector(
//             solution_x_.begin() + 11 * L,
//             solution_x_.begin() + 12 * L
//         );

//         std::vector<std::vector<double>> px_2d(1, px_vector);
//         std::vector<std::vector<double>> py_2d(1, py_vector);
//         std::vector<std::vector<double>> pz_2d(1, pz_vector);
//         std::vector<std::vector<double>> psi_2d(1, psi_vector);

//         std::vector<std::vector<double>> vx_2d(1, vx_vector);
//         std::vector<std::vector<double>> vy_2d(1, vy_vector);
//         std::vector<std::vector<double>> vz_2d(1, vz_vector);
//         std::vector<std::vector<double>> w_2d(1, w_vector);

//         std::vector<std::vector<double>> ax_2d(1, ax_vector);
//         std::vector<std::vector<double>> ay_2d(1, ay_vector);
//         std::vector<std::vector<double>> az_2d(1, az_vector);
//         std::vector<std::vector<double>> aw_2d(1, aw_vector);

//         auto px_real =
//             BernsteinPoly(px_2d, t_real, 0.0, tf_);

//         auto py_real =
//             BernsteinPoly(py_2d, t_real, 0.0, tf_);

//         auto pz_real =
//             BernsteinPoly(pz_2d, t_real, 0.0, tf_);

//         auto psi_real =
//             BernsteinPoly(psi_2d, t_real, 0.0, tf_);

//         auto vx_real =
//             BernsteinPoly(vx_2d, t_real, 0.0, tf_);

//         auto vy_real =
//             BernsteinPoly(vy_2d, t_real, 0.0, tf_);

//         auto vz_real =
//             BernsteinPoly(vz_2d, t_real, 0.0, tf_);

//         auto w_real =
//             BernsteinPoly(w_2d, t_real, 0.0, tf_);

//         auto ax_real =
//             BernsteinPoly(ax_2d, t_real, 0.0, tf_);

//         auto ay_real =
//             BernsteinPoly(ay_2d, t_real, 0.0, tf_);

//         auto az_real =
//             BernsteinPoly(az_2d, t_real, 0.0, tf_);

//         auto aw_real =
//             BernsteinPoly(aw_2d, t_real, 0.0, tf_);

//         auto flatten =
//             [](const std::vector<std::vector<double>>& input) {
//                 std::vector<double> output;

//                 for (const auto& row : input) {
//                     output.insert(output.end(), row.begin(), row.end());
//                 }

//                 return output;
//             };

//         writeToCSV(t_real, flatten(px_real), "px_real.csv");
//         writeToCSV(bebot_.getNodes(), px_vector, "px_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(py_real), "py_real.csv");
//         writeToCSV(bebot_.getNodes(), py_vector, "py_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(pz_real), "pz_real.csv");
//         writeToCSV(bebot_.getNodes(), pz_vector, "pz_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(psi_real), "psi_real.csv");
//         writeToCSV(bebot_.getNodes(), psi_vector, "psi_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(vx_real), "vx_real.csv");
//         writeToCSV(bebot_.getNodes(), vx_vector, "vx_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(vy_real), "vy_real.csv");
//         writeToCSV(bebot_.getNodes(), vy_vector, "vy_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(vz_real), "vz_real.csv");
//         writeToCSV(bebot_.getNodes(), vz_vector, "vz_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(w_real), "w_real.csv");
//         writeToCSV(bebot_.getNodes(), w_vector, "w_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(ax_real), "ax_real.csv");
//         writeToCSV(bebot_.getNodes(), ax_vector, "ax_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(ay_real), "ay_real.csv");
//         writeToCSV(bebot_.getNodes(), ay_vector, "ay_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(az_real), "az_real.csv");
//         writeToCSV(bebot_.getNodes(), az_vector, "az_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(aw_real), "aw_real.csv");
//         writeToCSV(bebot_.getNodes(), aw_vector, "aw_controlpoints_real.csv");
//     }

//     const std::vector<Number>& get_solution_x() const {
//         return solution_x_;
//     }

//     Number get_final_obj_value() const {
//         return final_obj_value_;
//     }

// private:
//     int N_;
//     double tf_;

//     double px_max_;
//     double px_min_;
//     double py_max_;
//     double py_min_;
//     double pz_max_;
//     double pz_min_;
//     double psi_max_;
//     double psi_min_;

//     double v_max_;
//     double w_max_;
//     double a_max_;
//     double aw_max_;

//     double px_cur_;
//     double py_cur_;
//     double pz_cur_;
//     double psi_cur_;

//     double vx_cur_;
//     double vy_cur_;
//     double vz_cur_;
//     double w_cur_;

//     double pxf_;
//     double pyf_;
//     double pzf_;
//     double psif_;

//     std::array<double, 1> obs_x_;
//     std::array<double, 1> obs_y_;
//     double obs_sep_;
//     int obs_deg_elev_extra_;

//     Bebot bebot_;

//     std::vector<double> tau_vec_;
//     std::vector<double> r_tau_;
//     std::vector<double> r_tau_row_;

//     std::vector<Number> solution_u_;
//     std::vector<Number> solution_x2_;
//     std::vector<Number> solution_x_;

//     Number final_obj_value_;

//     std::vector<double> final_time_;

//     std::vector<std::vector<double>> bernsteinpoly_resultu_;
//     std::vector<std::vector<double>> bernsteinpoly_resultx2_;
//     std::vector<std::vector<double>> bernsteinpoly_resultz_;

// public:
//     const std::vector<std::vector<double>>& get_bernsteinpoly_result() const {
//         return bernsteinpoly_resultz_;
//     }
// };

// extern "C" {
//     PointSetProblem* create_point_set_problem(
//         int N,
//         double tf,
//         double px_max,
//         double px_min,
//         double py_max,
//         double py_min,
//         double pz_max,
//         double pz_min,
//         double psi_max,
//         double psi_min,
//         double v_max,
//         double w_max,
//         double a_max,
//         double aw_max,
//         double px_cur,
//         double py_cur,
//         double pz_cur,
//         double psi_cur,
//         double vx_cur,
//         double vy_cur,
//         double vz_cur,
//         double w_cur,
//         double pxf,
//         double pyf,
//         double pzf,
//         double psif
//     ) {
//         return new PointSetProblem(
//             N,
//             tf,
//             px_max,
//             px_min,
//             py_max,
//             py_min,
//             pz_max,
//             pz_min,
//             psi_max,
//             psi_min,
//             v_max,
//             w_max,
//             a_max,
//             aw_max,
//             px_cur,
//             py_cur,
//             pz_cur,
//             psi_cur,
//             vx_cur,
//             vy_cur,
//             vz_cur,
//             w_cur,
//             pxf,
//             pyf,
//             pzf,
//             psif
//         );
//     }

//     void solve_point_set_problem(PointSetProblem* problem) {
//         SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

//         app->Options()->SetStringValue("linear_solver", "ma57");
//         app->Options()->SetStringValue("mu_strategy", "adaptive");

//         app->Options()->SetStringValue(
//             "gradient_approximation",
//             "finite-difference-values"
//         );

//         app->Options()->SetStringValue(
//             "jacobian_approximation",
//             "finite-difference-values"
//         );

//         app->Options()->SetStringValue(
//             "hessian_approximation",
//             "limited-memory"
//         );

//         app->Options()->SetIntegerValue("max_iter", 400);

//         app->Options()->SetNumericValue("tol", 1e-4);
//         app->Options()->SetNumericValue("constr_viol_tol", 1e-4);
//         app->Options()->SetNumericValue("obj_scaling_factor", 1e-4);

//         app->Options()->SetIntegerValue("print_level", 3);

//         app->RethrowNonIpoptException(true);

//         ApplicationReturnStatus status = app->Initialize();

//         if (status != Solve_Succeeded) {
//             std::cerr << "IPOPT initialization failed!" << std::endl;
//             return;
//         }

//         status = app->OptimizeTNLP(problem);

//         if (status == Solve_Succeeded || status == Solved_To_Acceptable_Level) {
//             std::cout << "Optimization succeeded!" << std::endl;

//             const auto& solution_x = problem->get_solution_x();

//             std::cout << "Optimal Solution (x): ";

//             for (Index i = 0; i < static_cast<Index>(solution_x.size()); i++) {
//                 std::cout << solution_x[i] << " ";
//             }

//             std::cout << std::endl;
//         } else {
//             std::cerr << "Optimization failed with status "
//                       << status << std::endl;
//         }
//     }

//     void get_solution(PointSetProblem* problem, double* solution, int n) {
//         const std::vector<double>& sol = problem->get_solution_x();

//         std::cout << "[DEBUG] get_solution called. Vector size: "
//                   << sol.size() << std::endl;

//         std::copy(sol.begin(), sol.end(), solution);
//     }

//     double get_final_objective_value(PointSetProblem* problem) {
//         return problem->get_final_obj_value();
//     }

//     void destroy_point_set_problem(PointSetProblem* problem) {
//         delete problem;
//     }
// }

////////////////////////////////////////////////////////////////////////////////////////////////////////
// two cylinder obstacles no memory issue


// #include "../../../../Ipopt_ma57_solver/src/Interfaces/IpIpoptApplication.hpp"
// #include "../../../../Ipopt_ma57_solver/src/Interfaces/IpTNLP.hpp"

// #include <cmath>
// #include <iostream>
// #include <fstream>
// #include <vector>
// #include <iomanip>
// #include <array>
// #include <string>
// #include <limits>
// #include <algorithm>
// #include <sstream>

// #include "../../../../include/bebot.h"
// #include "../../../../include/bernsteinpoly.h"
// #include "../../../../include/bernsteinproduct.h"
// #include "../../../../include/degelevmatrix.h"

// #include "mkl.h"

// using namespace Ipopt;

// class PointSetProblem : public Ipopt::TNLP {
// public:
//     PointSetProblem(
//         int N,
//         double tf,
//         double px_max,
//         double px_min,
//         double py_max,
//         double py_min,
//         double pz_max,
//         double pz_min,
//         double psi_max,
//         double psi_min,
//         double v_max,
//         double w_max,
//         double a_max,
//         double aw_max,
//         double px_cur,
//         double py_cur,
//         double pz_cur,
//         double psi_cur,
//         double vx_cur,
//         double vy_cur,
//         double vz_cur,
//         double w_cur,
//         double pxf,
//         double pyf,
//         double pzf,
//         double psif
//     )
//         : N_(N),
//           tf_(tf),
//           px_max_(px_max),
//           px_min_(px_min),
//           py_max_(py_max),
//           py_min_(py_min),
//           pz_max_(pz_max),
//           pz_min_(pz_min),
//           psi_max_(psi_max),
//           psi_min_(psi_min),
//           v_max_(v_max),
//           w_max_(w_max),
//           a_max_(a_max),
//           aw_max_(aw_max),
//           px_cur_(px_cur),
//           py_cur_(py_cur),
//           pz_cur_(pz_cur),
//           psi_cur_(psi_cur),
//           vx_cur_(vx_cur),
//           vy_cur_(vy_cur),
//           vz_cur_(vz_cur),
//           w_cur_(w_cur),
//           pxf_(pxf),
//           pyf_(pyf),
//           pzf_(pzf),
//           psif_(psif),
//           obs_x_{{-0.6, 0.6}},
//           obs_y_{{ 2.0, 1.8}},
//           obs_sep_(0.3),
//           obs_deg_elev_extra_(5),
//           bebot_(N, tf_) {
//         bebot_.calculate();

//         const int L = N_ + 1;

//         tau_vec_.resize(L);
//         r_tau_.resize(L);
//         r_tau_row_.resize(L);

//         if (L == 1) {
//             const double tau0 = -1.0;
//             const double denom = 1.0 - tau0;

//             tau_vec_[0] = tau0;
//             r_tau_[0] = 2.0 / (denom * denom);
//             r_tau_row_[0] = r_tau_[0];
//         } else {
//             const double tau0 = -1.0;
//             const double tau1 = 0.9;
//             const double step = (tau1 - tau0) / static_cast<double>(N_);

//             for (int i = 0; i < L; ++i) {
//                 const double tau = tau0 + step * static_cast<double>(i);
//                 const double denom = 1.0 - tau;

//                 tau_vec_[i] = tau;
//                 r_tau_[i] = 2.0 / (denom * denom);
//                 r_tau_row_[i] = r_tau_[i];
//             }
//         }
//     }

//     void writeToCSV(
//         const std::vector<double>& times,
//         const std::vector<double>& values,
//         const std::string& filename
//     ) {
//         std::ofstream outFile(filename);

//         if (!outFile.is_open()) {
//             std::cerr << "Failed to open file: " << filename << std::endl;
//             return;
//         }

//         outFile << "Time,Value\n";

//         for (size_t i = 0; i < times.size(); ++i) {
//             outFile << std::fixed << std::setprecision(6)
//                     << times[i] << "," << values[i] << "\n";
//         }

//         outFile.close();
//     }

//     virtual bool get_nlp_info(
//         Index& n,
//         Index& m,
//         Index& nnz_jac_g,
//         Index& nnz_h_lag,
//         IndexStyleEnum& index_style
//     ) {
//         const int L = N_ + 1;

//         n = 12 * L;

//         /*
//             Existing constraints:
//             8 * (N + 1): dynamics equality constraints
//             1 * (N + 1): speed inequality constraints
//             1 * (N + 1): acceleration inequality constraints

//             Obstacle constraints:
//             For each obstacle:
//                 dist2obs_square has degree 2N
//                 then degree-elevated to 2N + obs_deg_elev_extra_
//                 number of coefficients = 2N + obs_deg_elev_extra_ + 1

//             Total:
//                 m = 10 * (N + 1)
//                     + number_of_obstacles * (2N + obs_deg_elev_extra_ + 1)
//         */

//         const int obs_degree = 2 * N_ + obs_deg_elev_extra_;
//         const int obs_constraints =
//             static_cast<int>(obs_x_.size()) * (obs_degree + 1);

//         m = 10 * L + obs_constraints;

//         nnz_jac_g = n * m;
//         nnz_h_lag = 0;
//         index_style = TNLP::C_STYLE;

//         return true;
//     }

//     virtual bool get_bounds_info(
//         Index n,
//         Number* x_l,
//         Number* x_u,
//         Index m,
//         Number* g_l,
//         Number* g_u
//     ) {
//         const int L = N_ + 1;

//         std::vector<double> x_lower(
//             n,
//             -std::numeric_limits<double>::infinity()
//         );

//         std::vector<double> x_upper(
//             n,
//             std::numeric_limits<double>::infinity()
//         );

//         // px
//         for (int i = 1; i < L; ++i) {
//             x_lower[i] = px_min_;
//             x_upper[i] = px_max_;
//         }
//         x_lower[0] = px_cur_;
//         x_upper[0] = px_cur_;

//         // py
//         for (int i = L + 1; i < 2 * L; ++i) {
//             x_lower[i] = py_min_;
//             x_upper[i] = py_max_;
//         }
//         x_lower[L] = py_cur_;
//         x_upper[L] = py_cur_;

//         // pz
//         for (int i = 2 * L + 1; i < 3 * L; ++i) {
//             x_lower[i] = pz_min_;
//             x_upper[i] = pz_max_;
//         }
//         x_lower[2 * L] = pz_cur_;
//         x_upper[2 * L] = pz_cur_;

//         // psi
//         for (int i = 3 * L + 1; i < 4 * L; ++i) {
//             x_lower[i] = psi_min_;
//             x_upper[i] = psi_max_;
//         }
//         x_lower[3 * L] = psi_cur_;
//         x_upper[3 * L] = psi_cur_;

//         // vx
//         for (int i = 4 * L + 1; i < 5 * L; ++i) {
//             x_lower[i] = -v_max_;
//             x_upper[i] =  v_max_;
//         }
//         x_lower[4 * L] = vx_cur_;
//         x_upper[4 * L] = vx_cur_;

//         // vy
//         for (int i = 5 * L + 1; i < 6 * L; ++i) {
//             x_lower[i] = -v_max_;
//             x_upper[i] =  v_max_;
//         }
//         x_lower[5 * L] = vy_cur_;
//         x_upper[5 * L] = vy_cur_;

//         // vz
//         for (int i = 6 * L + 1; i < 7 * L; ++i) {
//             x_lower[i] = -v_max_;
//             x_upper[i] =  v_max_;
//         }
//         x_lower[6 * L] = vz_cur_;
//         x_upper[6 * L] = vz_cur_;

//         // w
//         for (int i = 7 * L + 1; i < 8 * L; ++i) {
//             x_lower[i] = -w_max_;
//             x_upper[i] =  w_max_;
//         }
//         x_lower[7 * L] = w_cur_;
//         x_upper[7 * L] = w_cur_;

//         // ax
//         for (int i = 8 * L; i < 9 * L; ++i) {
//             x_lower[i] = -a_max_;
//             x_upper[i] =  a_max_;
//         }

//         // ay
//         for (int i = 9 * L; i < 10 * L; ++i) {
//             x_lower[i] = -a_max_;
//             x_upper[i] =  a_max_;
//         }

//         // az
//         for (int i = 10 * L; i < 11 * L; ++i) {
//             x_lower[i] = -a_max_;
//             x_upper[i] =  a_max_;
//         }

//         // aw
//         for (int i = 11 * L; i < 12 * L; ++i) {
//             x_lower[i] = -aw_max_;
//             x_upper[i] =  aw_max_;
//         }

//         std::copy(x_lower.begin(), x_lower.end(), x_l);
//         std::copy(x_upper.begin(), x_upper.end(), x_u);

//         // Dynamics equality constraints
//         for (int i = 0; i < 8 * L; ++i) {
//             g_l[i] = 0.0;
//             g_u[i] = 0.0;
//         }

//         // Speed inequality constraints: c_speed <= 0
//         for (int i = 8 * L; i < 9 * L; ++i) {
//             g_l[i] = -std::numeric_limits<double>::infinity();
//             g_u[i] = 0.0;
//         }

//         // Acceleration inequality constraints: c_accel <= 0
//         for (int i = 9 * L; i < 10 * L; ++i) {
//             g_l[i] = -std::numeric_limits<double>::infinity();
//             g_u[i] = 0.0;
//         }

//         // Obstacle inequality constraints: c_obs <= 0
//         const int obs_start = 10 * L;

//         for (Index i = obs_start; i < m; ++i) {
//             g_l[i] = -std::numeric_limits<double>::infinity();
//             g_u[i] = 0.0;
//         }

//         return true;
//     }

//     virtual bool get_starting_point(
//         Index n,
//         bool init_x,
//         Number* x,
//         bool init_z,
//         Number* z_L,
//         Number* z_U,
//         Index m,
//         bool init_lambda,
//         Number* lambda
//     ) {
//         for (Index i = 0; i < n; ++i) {
//             x[i] = 1.0;
//         }

//         return true;
//     }

//     virtual bool eval_f(
//         Index n,
//         const Number* x,
//         bool new_x,
//         Number& obj_value
//     ) {
//         const int L = N_ + 1;

//         const double w_p   = 0.01;
//         const double w_psi = 1.0;
//         const double w_a   = 0.01;
//         const double w_aw  = 0.001;

//         obj_value = 0.0;

//         std::vector<double> pxf_vector(L, pxf_);
//         std::vector<double> pyf_vector(L, pyf_);
//         std::vector<double> pzf_vector(L, pzf_);
//         std::vector<double> psif_vector(L, psif_);

//         std::vector<double> px_vector(x, x + L);
//         std::vector<double> py_vector(x + 1 * L, x + 2 * L);
//         std::vector<double> pz_vector(x + 2 * L, x + 3 * L);
//         std::vector<double> psi_vector(x + 3 * L, x + 4 * L);

//         std::vector<double> ax_vector(x + 8  * L, x + 9  * L);
//         std::vector<double> ay_vector(x + 9  * L, x + 10 * L);
//         std::vector<double> az_vector(x + 10 * L, x + 11 * L);
//         std::vector<double> aw_vector(x + 11 * L, x + 12 * L);

//         std::vector<double> px_diff(L);
//         std::vector<double> py_diff(L);
//         std::vector<double> pz_diff(L);
//         std::vector<double> psi_diff(L);

//         vdSub(L, px_vector.data(),  pxf_vector.data(),  px_diff.data());
//         vdSub(L, py_vector.data(),  pyf_vector.data(),  py_diff.data());
//         vdSub(L, pz_vector.data(),  pzf_vector.data(),  pz_diff.data());
//         vdSub(L, psi_vector.data(), psif_vector.data(), psi_diff.data());

//         std::vector<double> px_diff_sqr(L);
//         std::vector<double> py_diff_sqr(L);
//         std::vector<double> pz_diff_sqr(L);
//         std::vector<double> psi_diff_sqr(L);

//         vdSqr(L, px_diff.data(),  px_diff_sqr.data());
//         vdSqr(L, py_diff.data(),  py_diff_sqr.data());
//         vdSqr(L, pz_diff.data(),  pz_diff_sqr.data());
//         vdSqr(L, psi_diff.data(), psi_diff_sqr.data());

//         std::vector<double> px_weighted(L);
//         std::vector<double> py_weighted(L);
//         std::vector<double> pz_weighted(L);
//         std::vector<double> psi_weighted(L);

//         vdMul(L, r_tau_.data(), px_diff_sqr.data(),  px_weighted.data());
//         vdMul(L, r_tau_.data(), py_diff_sqr.data(),  py_weighted.data());
//         vdMul(L, r_tau_.data(), pz_diff_sqr.data(),  pz_weighted.data());
//         vdMul(L, r_tau_.data(), psi_diff_sqr.data(), psi_weighted.data());

//         std::vector<double> ones(L, 1.0);

//         const double sum_px_state =
//             cblas_ddot(L, px_weighted.data(), 1, ones.data(), 1);

//         const double sum_py_state =
//             cblas_ddot(L, py_weighted.data(), 1, ones.data(), 1);

//         const double sum_pz_state =
//             cblas_ddot(L, pz_weighted.data(), 1, ones.data(), 1);

//         const double sum_psi_state =
//             cblas_ddot(L, psi_weighted.data(), 1, ones.data(), 1);

//         const double state_term =
//             w_p   * sum_px_state
//           + w_p   * sum_py_state
//           + w_p   * sum_pz_state
//           + w_psi * sum_psi_state;

//         std::vector<double> ax_sqr(L);
//         std::vector<double> ay_sqr(L);
//         std::vector<double> az_sqr(L);
//         std::vector<double> aw_sqr(L);

//         vdSqr(L, ax_vector.data(), ax_sqr.data());
//         vdSqr(L, ay_vector.data(), ay_sqr.data());
//         vdSqr(L, az_vector.data(), az_sqr.data());
//         vdSqr(L, aw_vector.data(), aw_sqr.data());

//         const double sum_ax =
//             cblas_ddot(L, ax_sqr.data(), 1, ones.data(), 1);

//         const double sum_ay =
//             cblas_ddot(L, ay_sqr.data(), 1, ones.data(), 1);

//         const double sum_az =
//             cblas_ddot(L, az_sqr.data(), 1, ones.data(), 1);

//         const double sum_aw =
//             cblas_ddot(L, aw_sqr.data(), 1, ones.data(), 1);

//         const double control_term =
//             w_a  * sum_ax
//           + w_a  * sum_ay
//           + w_a  * sum_az
//           + w_aw * sum_aw;

//         obj_value = state_term + control_term;

//         return true;
//     }

//     virtual bool eval_g(
//         Index n,
//         const Number* x,
//         bool new_x,
//         Index m,
//         Number* g
//     ) {
//         const int L = N_ + 1;

//         /*
//             Reuse already-computed BeBOT member.
//             bebot_.calculate() is called once in constructor.
//         */
//         const auto& Dm = bebot_.getDifferentiationMatrix();

//         std::vector<double> px_vector(x, x + L);
//         std::vector<double> py_vector(x + 1 * L, x + 2 * L);
//         std::vector<double> pz_vector(x + 2 * L, x + 3 * L);
//         std::vector<double> psi_vector(x + 3 * L, x + 4 * L);

//         std::vector<double> vx_vector(x + 4 * L, x + 5 * L);
//         std::vector<double> vy_vector(x + 5 * L, x + 6 * L);
//         std::vector<double> vz_vector(x + 6 * L, x + 7 * L);
//         std::vector<double> w_vector(x + 7 * L, x + 8 * L);

//         std::vector<double> ax_vector(x + 8  * L, x + 9  * L);
//         std::vector<double> ay_vector(x + 9  * L, x + 10 * L);
//         std::vector<double> az_vector(x + 10 * L, x + 11 * L);
//         std::vector<double> aw_vector(x + 11 * L, x + 12 * L);

//         std::vector<double> dyn1(L);
//         std::vector<double> dyn2(L);
//         std::vector<double> dyn3(L);
//         std::vector<double> dyn4(L);
//         std::vector<double> dyn5(L);
//         std::vector<double> dyn6(L);
//         std::vector<double> dyn7(L);
//         std::vector<double> dyn8(L);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, px_vector.data(), 1, 0.0, dyn1.data(), 1);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, py_vector.data(), 1, 0.0, dyn2.data(), 1);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, pz_vector.data(), 1, 0.0, dyn3.data(), 1);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, psi_vector.data(), 1, 0.0, dyn4.data(), 1);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, vx_vector.data(), 1, 0.0, dyn5.data(), 1);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, vy_vector.data(), 1, 0.0, dyn6.data(), 1);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, vz_vector.data(), 1, 0.0, dyn7.data(), 1);

//         cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
//                     Dm.data(), L, w_vector.data(), 1, 0.0, dyn8.data(), 1);

//         std::vector<double> px_rhs_scaled(L);
//         std::vector<double> py_rhs_scaled(L);
//         std::vector<double> pz_rhs_scaled(L);
//         std::vector<double> psi_rhs_scaled(L);

//         std::vector<double> vx_rhs_scaled(L);
//         std::vector<double> vy_rhs_scaled(L);
//         std::vector<double> vz_rhs_scaled(L);
//         std::vector<double> w_rhs_scaled(L);

//         vdMul(L, vx_vector.data(), r_tau_row_.data(), px_rhs_scaled.data());
//         vdMul(L, vy_vector.data(), r_tau_row_.data(), py_rhs_scaled.data());
//         vdMul(L, vz_vector.data(), r_tau_row_.data(), pz_rhs_scaled.data());
//         vdMul(L, w_vector.data(),  r_tau_row_.data(), psi_rhs_scaled.data());

//         vdMul(L, ax_vector.data(), r_tau_row_.data(), vx_rhs_scaled.data());
//         vdMul(L, ay_vector.data(), r_tau_row_.data(), vy_rhs_scaled.data());
//         vdMul(L, az_vector.data(), r_tau_row_.data(), vz_rhs_scaled.data());
//         vdMul(L, aw_vector.data(), r_tau_row_.data(), w_rhs_scaled.data());

//         std::vector<double> res_px(L);
//         std::vector<double> res_py(L);
//         std::vector<double> res_pz(L);
//         std::vector<double> res_psi(L);

//         std::vector<double> res_vx(L);
//         std::vector<double> res_vy(L);
//         std::vector<double> res_vz(L);
//         std::vector<double> res_w(L);

//         vdSub(L, dyn1.data(), px_rhs_scaled.data(),  res_px.data());
//         vdSub(L, dyn2.data(), py_rhs_scaled.data(),  res_py.data());
//         vdSub(L, dyn3.data(), pz_rhs_scaled.data(),  res_pz.data());
//         vdSub(L, dyn4.data(), psi_rhs_scaled.data(), res_psi.data());

//         vdSub(L, dyn5.data(), vx_rhs_scaled.data(), res_vx.data());
//         vdSub(L, dyn6.data(), vy_rhs_scaled.data(), res_vy.data());
//         vdSub(L, dyn7.data(), vz_rhs_scaled.data(), res_vz.data());
//         vdSub(L, dyn8.data(), w_rhs_scaled.data(),  res_w.data());

//         const double v_max2 = v_max_ * v_max_;
//         const double a_max2 = a_max_ * a_max_;

//         std::vector<double> vx2(L);
//         std::vector<double> vy2(L);
//         std::vector<double> vz2(L);

//         std::vector<double> ax2(L);
//         std::vector<double> ay2(L);
//         std::vector<double> az2(L);

//         vdSqr(L, vx_vector.data(), vx2.data());
//         vdSqr(L, vy_vector.data(), vy2.data());
//         vdSqr(L, vz_vector.data(), vz2.data());

//         vdSqr(L, ax_vector.data(), ax2.data());
//         vdSqr(L, ay_vector.data(), ay2.data());
//         vdSqr(L, az_vector.data(), az2.data());

//         std::vector<double> v2_sum(L);
//         std::vector<double> a2_sum(L);

//         vdAdd(L, vx2.data(), vy2.data(), v2_sum.data());
//         vdAdd(L, v2_sum.data(), vz2.data(), v2_sum.data());

//         vdAdd(L, ax2.data(), ay2.data(), a2_sum.data());
//         vdAdd(L, a2_sum.data(), az2.data(), a2_sum.data());

//         std::vector<double> v_max2_vec(L, v_max2);
//         std::vector<double> a_max2_vec(L, a_max2);

//         std::vector<double> c_speed(L);
//         std::vector<double> c_accel(L);

//         vdSub(L, v2_sum.data(), v_max2_vec.data(), c_speed.data());
//         vdSub(L, a2_sum.data(), a_max2_vec.data(), c_accel.data());

//         // ------------------------------------------------------------
//         // Obstacle avoidance constraints
//         //
//         // dist2obs_square =
//         //     BernsteinProduct(px - ox, px - ox)
//         //   + BernsteinProduct(py - oy, py - oy)
//         //
//         // c_obs = -dist2obs_square_elev + sep^2 <= 0
//         //
//         // Therefore:
//         //     dist2obs_square_elev >= sep^2
//         // ------------------------------------------------------------

//         auto add_vectors =
//             [](const std::vector<double>& a,
//             const std::vector<double>& b) -> std::vector<double> {
//                 const size_t size = std::min(a.size(), b.size());
//                 std::vector<double> out(size, 0.0);

//                 for (size_t i = 0; i < size; ++i) {
//                     out[i] = a[i] + b[i];
//                 }

//                 return out;
//             };

//         auto binom_ld =
//             [](int n, int k) -> long double {
//                 if (k < 0 || k > n) {
//                     return 0.0L;
//                 }

//                 if (k == 0 || k == n) {
//                     return 1.0L;
//                 }

//                 if (k > n - k) {
//                     k = n - k;
//                 }

//                 long double result = 1.0L;

//                 for (int i = 1; i <= k; ++i) {
//                     result *= static_cast<long double>(n - k + i);
//                     result /= static_cast<long double>(i);
//                 }

//                 return result;
//             };

//         auto degree_elevate =
//             [&](const std::vector<double>& cp,
//                 int target_degree) -> std::vector<double> {
//                 const int current_degree = static_cast<int>(cp.size()) - 1;

//                 if (target_degree <= current_degree) {
//                     return cp;
//                 }

//                 const int r = target_degree - current_degree;

//                 std::vector<double> elevated(target_degree + 1, 0.0);

//                 for (int j = 0; j <= target_degree; ++j) {
//                     const int i_min = std::max(0, j - r);
//                     const int i_max = std::min(current_degree, j);

//                     long double value = 0.0L;

//                     for (int i = i_min; i <= i_max; ++i) {
//                         const long double coeff =
//                             binom_ld(current_degree, i)
//                         * binom_ld(r, j - i)
//                         / binom_ld(target_degree, j);

//                         value += coeff * static_cast<long double>(cp[i]);
//                     }

//                     elevated[j] = static_cast<double>(value);
//                 }

//                 return elevated;
//             };

//         std::vector<double> c_obs_all;
//         const double sep2 = obs_sep_ * obs_sep_;

//         for (size_t obs_idx = 0; obs_idx < obs_x_.size(); ++obs_idx) {
//             const double ox = obs_x_[obs_idx];
//             const double oy = obs_y_[obs_idx];

//             std::vector<double> dx(L);
//             std::vector<double> dy(L);

//             for (int i = 0; i < L; ++i) {
//                 dx[i] = px_vector[i] - ox;
//                 dy[i] = py_vector[i] - oy;
//             }

//             std::vector<double> dx2 = BernsteinProduct(dx, dx);
//             std::vector<double> dy2 = BernsteinProduct(dy, dy);

//             std::vector<double> dist2obs_square =
//                 add_vectors(dx2, dy2);

//             const int deg_current =
//                 static_cast<int>(dist2obs_square.size()) - 1;

//             const int deg_target =
//                 deg_current + obs_deg_elev_extra_;

//             std::vector<double> dist2obs_square_elev =
//                 degree_elevate(dist2obs_square, deg_target);

//             for (double coeff : dist2obs_square_elev) {
//                 c_obs_all.push_back(-coeff + sep2);
//             }
//         }

//         for (Index i = 0; i < L; ++i) {
//             g[0 * L + i] = res_px[i];
//             g[1 * L + i] = res_py[i];
//             g[2 * L + i] = res_pz[i];
//             g[3 * L + i] = res_psi[i];

//             g[4 * L + i] = res_vx[i];
//             g[5 * L + i] = res_vy[i];
//             g[6 * L + i] = res_vz[i];
//             g[7 * L + i] = res_w[i];

//             g[8 * L + i] = c_speed[i];
//             g[9 * L + i] = c_accel[i];
//         }

//         const Index obs_start = 10 * L;

//         for (Index i = 0; i < static_cast<Index>(c_obs_all.size()); ++i) {
//             g[obs_start + i] = c_obs_all[i];
//         }

//         return true;
//     }

//     virtual bool eval_jac_g(
//         Index n,
//         const Number* x,
//         bool new_x,
//         Index m,
//         Index nele_jac,
//         Index* iRow,
//         Index* jCol,
//         Number* values
//     ) {
//         if (values == NULL) {
//             for (Index i = 0; i < m; i++) {
//                 for (Index j = 0; j < n; j++) {
//                     iRow[i * n + j] = i;
//                     jCol[i * n + j] = j;
//                 }
//             }
//         }

//         return true;
//     }

//     virtual bool eval_grad_f(
//         Index n,
//         const Number* x,
//         bool new_x,
//         Number* grad_f
//     ) {
//         return true;
//     }

//     virtual void finalize_solution(
//         SolverReturn status,
//         Index n,
//         const Number* x,
//         const Number* z_L,
//         const Number* z_U,
//         Index m,
//         const Number* g,
//         const Number* lambda,
//         Number obj_value,
//         const IpoptData* ip_data,
//         IpoptCalculatedQuantities* ip_cq
//     ) {
//         const int L = N_ + 1;

//         solution_x_.resize(12 * L);

//         for (Index i = 0; i < 12 * L; ++i) {
//             solution_x_[i] = x[i];
//         }

//         final_obj_value_ = obj_value;

//         const int K = 1000;

//         std::vector<double> t_norm(K);
//         std::vector<double> t_real(K);

//         for (int i = 0; i < K; ++i) {
//             const double s =
//                 static_cast<double>(i) / static_cast<double>(K - 1);

//             t_norm[i] = s;
//             t_real[i] = s * tf_;
//         }

//         std::vector<double> px_vector(
//             solution_x_.begin(),
//             solution_x_.begin() + L
//         );

//         std::vector<double> py_vector(
//             solution_x_.begin() + 1 * L,
//             solution_x_.begin() + 2 * L
//         );

//         std::vector<double> pz_vector(
//             solution_x_.begin() + 2 * L,
//             solution_x_.begin() + 3 * L
//         );

//         std::vector<double> psi_vector(
//             solution_x_.begin() + 3 * L,
//             solution_x_.begin() + 4 * L
//         );

//         std::vector<double> vx_vector(
//             solution_x_.begin() + 4 * L,
//             solution_x_.begin() + 5 * L
//         );

//         std::vector<double> vy_vector(
//             solution_x_.begin() + 5 * L,
//             solution_x_.begin() + 6 * L
//         );

//         std::vector<double> vz_vector(
//             solution_x_.begin() + 6 * L,
//             solution_x_.begin() + 7 * L
//         );

//         std::vector<double> w_vector(
//             solution_x_.begin() + 7 * L,
//             solution_x_.begin() + 8 * L
//         );

//         std::vector<double> ax_vector(
//             solution_x_.begin() + 8 * L,
//             solution_x_.begin() + 9 * L
//         );

//         std::vector<double> ay_vector(
//             solution_x_.begin() + 9 * L,
//             solution_x_.begin() + 10 * L
//         );

//         std::vector<double> az_vector(
//             solution_x_.begin() + 10 * L,
//             solution_x_.begin() + 11 * L
//         );

//         std::vector<double> aw_vector(
//             solution_x_.begin() + 11 * L,
//             solution_x_.begin() + 12 * L
//         );

//         std::vector<std::vector<double>> px_2d(1, px_vector);
//         std::vector<std::vector<double>> py_2d(1, py_vector);
//         std::vector<std::vector<double>> pz_2d(1, pz_vector);
//         std::vector<std::vector<double>> psi_2d(1, psi_vector);

//         std::vector<std::vector<double>> vx_2d(1, vx_vector);
//         std::vector<std::vector<double>> vy_2d(1, vy_vector);
//         std::vector<std::vector<double>> vz_2d(1, vz_vector);
//         std::vector<std::vector<double>> w_2d(1, w_vector);

//         std::vector<std::vector<double>> ax_2d(1, ax_vector);
//         std::vector<std::vector<double>> ay_2d(1, ay_vector);
//         std::vector<std::vector<double>> az_2d(1, az_vector);
//         std::vector<std::vector<double>> aw_2d(1, aw_vector);

//         auto px_real =
//             BernsteinPoly(px_2d, t_real, 0.0, tf_);

//         auto py_real =
//             BernsteinPoly(py_2d, t_real, 0.0, tf_);

//         auto pz_real =
//             BernsteinPoly(pz_2d, t_real, 0.0, tf_);

//         auto psi_real =
//             BernsteinPoly(psi_2d, t_real, 0.0, tf_);

//         auto vx_real =
//             BernsteinPoly(vx_2d, t_real, 0.0, tf_);

//         auto vy_real =
//             BernsteinPoly(vy_2d, t_real, 0.0, tf_);

//         auto vz_real =
//             BernsteinPoly(vz_2d, t_real, 0.0, tf_);

//         auto w_real =
//             BernsteinPoly(w_2d, t_real, 0.0, tf_);

//         auto ax_real =
//             BernsteinPoly(ax_2d, t_real, 0.0, tf_);

//         auto ay_real =
//             BernsteinPoly(ay_2d, t_real, 0.0, tf_);

//         auto az_real =
//             BernsteinPoly(az_2d, t_real, 0.0, tf_);

//         auto aw_real =
//             BernsteinPoly(aw_2d, t_real, 0.0, tf_);

//         auto flatten =
//             [](const std::vector<std::vector<double>>& input) {
//                 std::vector<double> output;

//                 for (const auto& row : input) {
//                     output.insert(output.end(), row.begin(), row.end());
//                 }

//                 return output;
//             };

//         writeToCSV(t_real, flatten(px_real), "px_real.csv");
//         writeToCSV(bebot_.getNodes(), px_vector, "px_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(py_real), "py_real.csv");
//         writeToCSV(bebot_.getNodes(), py_vector, "py_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(pz_real), "pz_real.csv");
//         writeToCSV(bebot_.getNodes(), pz_vector, "pz_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(psi_real), "psi_real.csv");
//         writeToCSV(bebot_.getNodes(), psi_vector, "psi_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(vx_real), "vx_real.csv");
//         writeToCSV(bebot_.getNodes(), vx_vector, "vx_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(vy_real), "vy_real.csv");
//         writeToCSV(bebot_.getNodes(), vy_vector, "vy_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(vz_real), "vz_real.csv");
//         writeToCSV(bebot_.getNodes(), vz_vector, "vz_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(w_real), "w_real.csv");
//         writeToCSV(bebot_.getNodes(), w_vector, "w_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(ax_real), "ax_real.csv");
//         writeToCSV(bebot_.getNodes(), ax_vector, "ax_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(ay_real), "ay_real.csv");
//         writeToCSV(bebot_.getNodes(), ay_vector, "ay_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(az_real), "az_real.csv");
//         writeToCSV(bebot_.getNodes(), az_vector, "az_controlpoints_real.csv");

//         writeToCSV(t_real, flatten(aw_real), "aw_real.csv");
//         writeToCSV(bebot_.getNodes(), aw_vector, "aw_controlpoints_real.csv");
//     }

//     const std::vector<Number>& get_solution_x() const {
//         return solution_x_;
//     }

//     Number get_final_obj_value() const {
//         return final_obj_value_;
//     }

// private:
//     int N_;
//     double tf_;

//     double px_max_;
//     double px_min_;
//     double py_max_;
//     double py_min_;
//     double pz_max_;
//     double pz_min_;
//     double psi_max_;
//     double psi_min_;

//     double v_max_;
//     double w_max_;
//     double a_max_;
//     double aw_max_;

//     double px_cur_;
//     double py_cur_;
//     double pz_cur_;
//     double psi_cur_;

//     double vx_cur_;
//     double vy_cur_;
//     double vz_cur_;
//     double w_cur_;

//     double pxf_;
//     double pyf_;
//     double pzf_;
//     double psif_;

//     std::array<double, 2> obs_x_;
//     std::array<double, 2> obs_y_;
//     double obs_sep_;
//     int obs_deg_elev_extra_;

//     Bebot bebot_;

//     std::vector<double> tau_vec_;
//     std::vector<double> r_tau_;
//     std::vector<double> r_tau_row_;

//     std::vector<Number> solution_u_;
//     std::vector<Number> solution_x2_;
//     std::vector<Number> solution_x_;

//     Number final_obj_value_;

//     std::vector<double> final_time_;

//     std::vector<std::vector<double>> bernsteinpoly_resultu_;
//     std::vector<std::vector<double>> bernsteinpoly_resultx2_;
//     std::vector<std::vector<double>> bernsteinpoly_resultz_;

// public:
//     const std::vector<std::vector<double>>& get_bernsteinpoly_result() const {
//         return bernsteinpoly_resultz_;
//     }
// };

// extern "C" {
//     PointSetProblem* create_point_set_problem(
//         int N,
//         double tf,
//         double px_max,
//         double px_min,
//         double py_max,
//         double py_min,
//         double pz_max,
//         double pz_min,
//         double psi_max,
//         double psi_min,
//         double v_max,
//         double w_max,
//         double a_max,
//         double aw_max,
//         double px_cur,
//         double py_cur,
//         double pz_cur,
//         double psi_cur,
//         double vx_cur,
//         double vy_cur,
//         double vz_cur,
//         double w_cur,
//         double pxf,
//         double pyf,
//         double pzf,
//         double psif
//     ) {
//         PointSetProblem* problem = new PointSetProblem(
//             N,
//             tf,
//             px_max,
//             px_min,
//             py_max,
//             py_min,
//             pz_max,
//             pz_min,
//             psi_max,
//             psi_min,
//             v_max,
//             w_max,
//             a_max,
//             aw_max,
//             px_cur,
//             py_cur,
//             pz_cur,
//             psi_cur,
//             vx_cur,
//             vy_cur,
//             vz_cur,
//             w_cur,
//             pxf,
//             pyf,
//             pzf,
//             psif
//         );

//         /*
//             IMPORTANT LIFETIME FIX:

//             IPOPT uses SmartPtr reference counting for TNLP objects.
//             The external C caller receives a raw pointer, so we keep one
//             explicit reference for that external owner.

//             solve_point_set_problem() will create a temporary SmartPtr<TNLP>,
//             which adds/releases its own reference. This external reference
//             keeps the object alive after OptimizeTNLP() returns.

//             destroy_point_set_problem() releases this external reference.
//         */
//         problem->AddRef(nullptr);

//         return problem;
//     }

//     void solve_point_set_problem(PointSetProblem* problem) {
//         if (!problem) {
//             std::cerr << "solve_point_set_problem() received nullptr"
//                       << std::endl;
//             return;
//         }

//         SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

//         app->Options()->SetStringValue("linear_solver", "ma57");
//         app->Options()->SetStringValue("mu_strategy", "adaptive");

//         app->Options()->SetStringValue(
//             "gradient_approximation",
//             "finite-difference-values"
//         );

//         app->Options()->SetStringValue(
//             "jacobian_approximation",
//             "finite-difference-values"
//         );

//         app->Options()->SetStringValue(
//             "hessian_approximation",
//             "limited-memory"
//         );

//         app->Options()->SetIntegerValue("max_iter", 400);

//         app->Options()->SetNumericValue("tol", 1e-3);
//         app->Options()->SetNumericValue("constr_viol_tol", 1e-3);
//         app->Options()->SetNumericValue("obj_scaling_factor", 1e-3);

//         app->Options()->SetIntegerValue("print_level", 3);

//         app->RethrowNonIpoptException(true);

//         ApplicationReturnStatus status = app->Initialize();

//         if (status != Solve_Succeeded) {
//             std::cerr << "IPOPT initialization failed!" << std::endl;
//             return;
//         }

//         /*
//             IMPORTANT LIFETIME FIX:

//             Do not pass the raw pointer directly as OptimizeTNLP(problem).
//             Wrap it explicitly in SmartPtr<TNLP>. Because create_point_set_problem()
//             already called AddRef(), the object remains valid after this function
//             returns.
//         */
//         SmartPtr<TNLP> tnlp = problem;
//         status = app->OptimizeTNLP(tnlp);

//         if (status == Solve_Succeeded || status == Solved_To_Acceptable_Level) {
//             std::cout << "Optimization succeeded!" << std::endl;

//             const auto& solution_x = problem->get_solution_x();

//             std::cout << "Optimal Solution (x): ";

//             for (Index i = 0; i < static_cast<Index>(solution_x.size()); i++) {
//                 std::cout << solution_x[i] << " ";
//             }

//             std::cout << std::endl;
//         } else {
//             std::cerr << "Optimization failed with status "
//                       << status << std::endl;
//         }
//     }

//     void get_solution(PointSetProblem* problem, double* solution, int n) {
//         if (!problem) {
//             std::cerr << "get_solution() received nullptr problem"
//                       << std::endl;
//             return;
//         }

//         if (!solution) {
//             std::cerr << "get_solution() received nullptr solution buffer"
//                       << std::endl;
//             return;
//         }

//         if (n <= 0) {
//             std::cerr << "get_solution() received non-positive n = "
//                       << n << std::endl;
//             return;
//         }

//         const std::vector<double>& sol = problem->get_solution_x();

//         std::cout << "[DEBUG] get_solution called. Vector size: "
//                   << sol.size() << ", requested n: " << n << std::endl;

//         const int copy_n =
//             std::min(n, static_cast<int>(sol.size()));

//         std::copy(sol.begin(), sol.begin() + copy_n, solution);

//         if (copy_n < n) {
//             std::fill(solution + copy_n, solution + n, 0.0);

//             std::cerr << "[WARN] get_solution requested " << n
//                       << " values, but solution has only " << sol.size()
//                       << ". Remaining output entries filled with 0."
//                       << std::endl;
//         }
//     }

//     double get_final_objective_value(PointSetProblem* problem) {
//         if (!problem) {
//             std::cerr << "get_final_objective_value() received nullptr"
//                       << std::endl;
//             return std::numeric_limits<double>::quiet_NaN();
//         }

//         return problem->get_final_obj_value();
//     }

//     void destroy_point_set_problem(PointSetProblem* problem) {
//         if (problem) {
//             /*
//                 IMPORTANT LIFETIME FIX:

//                 Do not call delete problem.
//                 Release the external reference created in create_point_set_problem().
//             */
//             problem->ReleaseRef(nullptr);
//         }
//     }
// }

////////////////////////////////////////////////////////////////////////////////////////////////////////
// one cylinder obstacle, one shpere no memory issue

#include "../../../../Ipopt_ma57_solver/src/Interfaces/IpIpoptApplication.hpp"
#include "../../../../Ipopt_ma57_solver/src/Interfaces/IpTNLP.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <array>
#include <string>
#include <limits>
#include <algorithm>
#include <sstream>

#include "../../../../include/bebot.h"
#include "../../../../include/bernsteinpoly.h"
#include "../../../../include/bernsteinproduct.h"
#include "../../../../include/degelevmatrix.h"

#include "mkl.h"

using namespace Ipopt;

class PointSetProblem : public Ipopt::TNLP {
public:
    PointSetProblem(
        int N,
        double tf,
        double px_max,
        double px_min,
        double py_max,
        double py_min,
        double pz_max,
        double pz_min,
        double psi_max,
        double psi_min,
        double v_max,
        double w_max,
        double a_max,
        double aw_max,
        double px_cur,
        double py_cur,
        double pz_cur,
        double psi_cur,
        double vx_cur,
        double vy_cur,
        double vz_cur,
        double w_cur,
        double pxf,
        double pyf,
        double pzf,
        double psif,

        // Cylinder obstacle:
        // vertical cylinder-like obstacle, same style as before:
        // (px - cyl_x)^2 + (py - cyl_y)^2 >= cyl_radius^2
        double cyl_x,
        double cyl_y,
        double cyl_radius,

        // Sphere obstacle:
        // (px - sphere_x)^2 + (py - sphere_y)^2 + (pz - sphere_z)^2 >= sphere_radius^2
        double sphere_x,
        double sphere_y,
        double sphere_z,
        double sphere_radius
    )
        : N_(N),
          tf_(tf),
          px_max_(px_max),
          px_min_(px_min),
          py_max_(py_max),
          py_min_(py_min),
          pz_max_(pz_max),
          pz_min_(pz_min),
          psi_max_(psi_max),
          psi_min_(psi_min),
          v_max_(v_max),
          w_max_(w_max),
          a_max_(a_max),
          aw_max_(aw_max),
          px_cur_(px_cur),
          py_cur_(py_cur),
          pz_cur_(pz_cur),
          psi_cur_(psi_cur),
          vx_cur_(vx_cur),
          vy_cur_(vy_cur),
          vz_cur_(vz_cur),
          w_cur_(w_cur),
          pxf_(pxf),
          pyf_(pyf),
          pzf_(pzf),
          psif_(psif),
          cyl_x_(cyl_x),
          cyl_y_(cyl_y),
          cyl_radius_(cyl_radius),
          sphere_x_(sphere_x),
          sphere_y_(sphere_y),
          sphere_z_(sphere_z),
          sphere_radius_(sphere_radius),
          obs_deg_elev_extra_(5),
          bebot_(N, tf_) {
        bebot_.calculate();

        const int L = N_ + 1;

        tau_vec_.resize(L);
        r_tau_.resize(L);
        r_tau_row_.resize(L);

        if (L == 1) {
            const double tau0 = -1.0;
            const double denom = 1.0 - tau0;

            tau_vec_[0] = tau0;
            r_tau_[0] = 2.0 / (denom * denom);
            r_tau_row_[0] = r_tau_[0];
        } else {
            const double tau0 = -1.0;
            const double tau1 = 0.9;
            const double step = (tau1 - tau0) / static_cast<double>(N_);

            for (int i = 0; i < L; ++i) {
                const double tau = tau0 + step * static_cast<double>(i);
                const double denom = 1.0 - tau;

                tau_vec_[i] = tau;
                r_tau_[i] = 2.0 / (denom * denom);
                r_tau_row_[i] = r_tau_[i];
            }
        }
    }

    void writeToCSV(
        const std::vector<double>& times,
        const std::vector<double>& values,
        const std::string& filename
    ) {
        std::ofstream outFile(filename);

        if (!outFile.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        outFile << "Time,Value\n";

        for (size_t i = 0; i < times.size(); ++i) {
            outFile << std::fixed << std::setprecision(6)
                    << times[i] << "," << values[i] << "\n";
        }

        outFile.close();
    }

    virtual bool get_nlp_info(
        Index& n,
        Index& m,
        Index& nnz_jac_g,
        Index& nnz_h_lag,
        IndexStyleEnum& index_style
    ) {
        const int L = N_ + 1;

        n = 12 * L;

        /*
            Existing constraints:
            8 * (N + 1): dynamics equality constraints
            1 * (N + 1): speed inequality constraints
            1 * (N + 1): acceleration inequality constraints

            Obstacle constraints:
            1 cylinder obstacle:
                degree 2N, elevated to 2N + obs_deg_elev_extra_
                number of coefficients = 2N + obs_deg_elev_extra_ + 1

            1 sphere obstacle:
                degree 2N, elevated to 2N + obs_deg_elev_extra_
                number of coefficients = 2N + obs_deg_elev_extra_ + 1

            Total:
                m = 10 * (N + 1)
                    + 2 * (2N + obs_deg_elev_extra_ + 1)
        */

        const int obs_degree = 2 * N_ + obs_deg_elev_extra_;
        const int obs_constraints = 2 * (obs_degree + 1);

        m = 10 * L + obs_constraints;

        nnz_jac_g = n * m;
        nnz_h_lag = 0;
        index_style = TNLP::C_STYLE;

        return true;
    }

    virtual bool get_bounds_info(
        Index n,
        Number* x_l,
        Number* x_u,
        Index m,
        Number* g_l,
        Number* g_u
    ) {
        const int L = N_ + 1;

        std::vector<double> x_lower(
            n,
            -std::numeric_limits<double>::infinity()
        );

        std::vector<double> x_upper(
            n,
            std::numeric_limits<double>::infinity()
        );

        // px
        for (int i = 1; i < L; ++i) {
            x_lower[i] = px_min_;
            x_upper[i] = px_max_;
        }
        x_lower[0] = px_cur_;
        x_upper[0] = px_cur_;

        // py
        for (int i = L + 1; i < 2 * L; ++i) {
            x_lower[i] = py_min_;
            x_upper[i] = py_max_;
        }
        x_lower[L] = py_cur_;
        x_upper[L] = py_cur_;

        // pz
        for (int i = 2 * L + 1; i < 3 * L; ++i) {
            x_lower[i] = pz_min_;
            x_upper[i] = pz_max_;
        }
        x_lower[2 * L] = pz_cur_;
        x_upper[2 * L] = pz_cur_;

        // psi
        for (int i = 3 * L + 1; i < 4 * L; ++i) {
            x_lower[i] = psi_min_;
            x_upper[i] = psi_max_;
        }
        x_lower[3 * L] = psi_cur_;
        x_upper[3 * L] = psi_cur_;

        // vx
        for (int i = 4 * L + 1; i < 5 * L; ++i) {
            x_lower[i] = -v_max_;
            x_upper[i] =  v_max_;
        }
        x_lower[4 * L] = vx_cur_;
        x_upper[4 * L] = vx_cur_;

        // vy
        for (int i = 5 * L + 1; i < 6 * L; ++i) {
            x_lower[i] = -v_max_;
            x_upper[i] =  v_max_;
        }
        x_lower[5 * L] = vy_cur_;
        x_upper[5 * L] = vy_cur_;

        // vz
        for (int i = 6 * L + 1; i < 7 * L; ++i) {
            x_lower[i] = -v_max_;
            x_upper[i] =  v_max_;
        }
        x_lower[6 * L] = vz_cur_;
        x_upper[6 * L] = vz_cur_;

        // w
        for (int i = 7 * L + 1; i < 8 * L; ++i) {
            x_lower[i] = -w_max_;
            x_upper[i] =  w_max_;
        }
        x_lower[7 * L] = w_cur_;
        x_upper[7 * L] = w_cur_;

        // ax
        for (int i = 8 * L; i < 9 * L; ++i) {
            x_lower[i] = -a_max_;
            x_upper[i] =  a_max_;
        }

        // ay
        for (int i = 9 * L; i < 10 * L; ++i) {
            x_lower[i] = -a_max_;
            x_upper[i] =  a_max_;
        }

        // az
        for (int i = 10 * L; i < 11 * L; ++i) {
            x_lower[i] = -a_max_;
            x_upper[i] =  a_max_;
        }

        // aw
        for (int i = 11 * L; i < 12 * L; ++i) {
            x_lower[i] = -aw_max_;
            x_upper[i] =  aw_max_;
        }

        std::copy(x_lower.begin(), x_lower.end(), x_l);
        std::copy(x_upper.begin(), x_upper.end(), x_u);

        // Dynamics equality constraints
        for (int i = 0; i < 8 * L; ++i) {
            g_l[i] = 0.0;
            g_u[i] = 0.0;
        }

        // Speed inequality constraints: c_speed <= 0
        for (int i = 8 * L; i < 9 * L; ++i) {
            g_l[i] = -std::numeric_limits<double>::infinity();
            g_u[i] = 0.0;
        }

        // Acceleration inequality constraints: c_accel <= 0
        for (int i = 9 * L; i < 10 * L; ++i) {
            g_l[i] = -std::numeric_limits<double>::infinity();
            g_u[i] = 0.0;
        }

        // Obstacle inequality constraints: c_obs <= 0
        const int obs_start = 10 * L;

        for (Index i = obs_start; i < m; ++i) {
            g_l[i] = -std::numeric_limits<double>::infinity();
            g_u[i] = 0.0;
        }

        return true;
    }

    virtual bool get_starting_point(
        Index n,
        bool init_x,
        Number* x,
        bool init_z,
        Number* z_L,
        Number* z_U,
        Index m,
        bool init_lambda,
        Number* lambda
    ) {
        for (Index i = 0; i < n; ++i) {
            x[i] = 1.0;
        }

        return true;
    }

    virtual bool eval_f(
        Index n,
        const Number* x,
        bool new_x,
        Number& obj_value
    ) {
        const int L = N_ + 1;

        const double w_p   = 0.01;
        const double w_psi = 1.0;
        const double w_a   = 0.01;
        const double w_aw  = 0.001;

        obj_value = 0.0;

        std::vector<double> pxf_vector(L, pxf_);
        std::vector<double> pyf_vector(L, pyf_);
        std::vector<double> pzf_vector(L, pzf_);
        std::vector<double> psif_vector(L, psif_);

        std::vector<double> px_vector(x, x + L);
        std::vector<double> py_vector(x + 1 * L, x + 2 * L);
        std::vector<double> pz_vector(x + 2 * L, x + 3 * L);
        std::vector<double> psi_vector(x + 3 * L, x + 4 * L);

        std::vector<double> ax_vector(x + 8  * L, x + 9  * L);
        std::vector<double> ay_vector(x + 9  * L, x + 10 * L);
        std::vector<double> az_vector(x + 10 * L, x + 11 * L);
        std::vector<double> aw_vector(x + 11 * L, x + 12 * L);

        std::vector<double> px_diff(L);
        std::vector<double> py_diff(L);
        std::vector<double> pz_diff(L);
        std::vector<double> psi_diff(L);

        vdSub(L, px_vector.data(),  pxf_vector.data(),  px_diff.data());
        vdSub(L, py_vector.data(),  pyf_vector.data(),  py_diff.data());
        vdSub(L, pz_vector.data(),  pzf_vector.data(),  pz_diff.data());
        vdSub(L, psi_vector.data(), psif_vector.data(), psi_diff.data());

        std::vector<double> px_diff_sqr(L);
        std::vector<double> py_diff_sqr(L);
        std::vector<double> pz_diff_sqr(L);
        std::vector<double> psi_diff_sqr(L);

        vdSqr(L, px_diff.data(),  px_diff_sqr.data());
        vdSqr(L, py_diff.data(),  py_diff_sqr.data());
        vdSqr(L, pz_diff.data(),  pz_diff_sqr.data());
        vdSqr(L, psi_diff.data(), psi_diff_sqr.data());

        std::vector<double> px_weighted(L);
        std::vector<double> py_weighted(L);
        std::vector<double> pz_weighted(L);
        std::vector<double> psi_weighted(L);

        vdMul(L, r_tau_.data(), px_diff_sqr.data(),  px_weighted.data());
        vdMul(L, r_tau_.data(), py_diff_sqr.data(),  py_weighted.data());
        vdMul(L, r_tau_.data(), pz_diff_sqr.data(),  pz_weighted.data());
        vdMul(L, r_tau_.data(), psi_diff_sqr.data(), psi_weighted.data());

        std::vector<double> ones(L, 1.0);

        const double sum_px_state =
            cblas_ddot(L, px_weighted.data(), 1, ones.data(), 1);

        const double sum_py_state =
            cblas_ddot(L, py_weighted.data(), 1, ones.data(), 1);

        const double sum_pz_state =
            cblas_ddot(L, pz_weighted.data(), 1, ones.data(), 1);

        const double sum_psi_state =
            cblas_ddot(L, psi_weighted.data(), 1, ones.data(), 1);

        const double state_term =
            w_p   * sum_px_state
          + w_p   * sum_py_state
          + w_p   * sum_pz_state
          + w_psi * sum_psi_state;

        std::vector<double> ax_sqr(L);
        std::vector<double> ay_sqr(L);
        std::vector<double> az_sqr(L);
        std::vector<double> aw_sqr(L);

        vdSqr(L, ax_vector.data(), ax_sqr.data());
        vdSqr(L, ay_vector.data(), ay_sqr.data());
        vdSqr(L, az_vector.data(), az_sqr.data());
        vdSqr(L, aw_vector.data(), aw_sqr.data());

        const double sum_ax =
            cblas_ddot(L, ax_sqr.data(), 1, ones.data(), 1);

        const double sum_ay =
            cblas_ddot(L, ay_sqr.data(), 1, ones.data(), 1);

        const double sum_az =
            cblas_ddot(L, az_sqr.data(), 1, ones.data(), 1);

        const double sum_aw =
            cblas_ddot(L, aw_sqr.data(), 1, ones.data(), 1);

        const double control_term =
            w_a  * sum_ax
          + w_a  * sum_ay
          + w_a  * sum_az
          + w_aw * sum_aw;

        obj_value = state_term + control_term;

        return true;
    }

    virtual bool eval_g(
        Index n,
        const Number* x,
        bool new_x,
        Index m,
        Number* g
    ) {
        const int L = N_ + 1;

        const auto& Dm = bebot_.getDifferentiationMatrix();

        std::vector<double> px_vector(x, x + L);
        std::vector<double> py_vector(x + 1 * L, x + 2 * L);
        std::vector<double> pz_vector(x + 2 * L, x + 3 * L);
        std::vector<double> psi_vector(x + 3 * L, x + 4 * L);

        std::vector<double> vx_vector(x + 4 * L, x + 5 * L);
        std::vector<double> vy_vector(x + 5 * L, x + 6 * L);
        std::vector<double> vz_vector(x + 6 * L, x + 7 * L);
        std::vector<double> w_vector(x + 7 * L, x + 8 * L);

        std::vector<double> ax_vector(x + 8  * L, x + 9  * L);
        std::vector<double> ay_vector(x + 9  * L, x + 10 * L);
        std::vector<double> az_vector(x + 10 * L, x + 11 * L);
        std::vector<double> aw_vector(x + 11 * L, x + 12 * L);

        std::vector<double> dyn1(L);
        std::vector<double> dyn2(L);
        std::vector<double> dyn3(L);
        std::vector<double> dyn4(L);
        std::vector<double> dyn5(L);
        std::vector<double> dyn6(L);
        std::vector<double> dyn7(L);
        std::vector<double> dyn8(L);

        cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
                    Dm.data(), L, px_vector.data(), 1, 0.0, dyn1.data(), 1);

        cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
                    Dm.data(), L, py_vector.data(), 1, 0.0, dyn2.data(), 1);

        cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
                    Dm.data(), L, pz_vector.data(), 1, 0.0, dyn3.data(), 1);

        cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
                    Dm.data(), L, psi_vector.data(), 1, 0.0, dyn4.data(), 1);

        cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
                    Dm.data(), L, vx_vector.data(), 1, 0.0, dyn5.data(), 1);

        cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
                    Dm.data(), L, vy_vector.data(), 1, 0.0, dyn6.data(), 1);

        cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
                    Dm.data(), L, vz_vector.data(), 1, 0.0, dyn7.data(), 1);

        cblas_dgemv(CblasColMajor, CblasTrans, L, L, 1.0,
                    Dm.data(), L, w_vector.data(), 1, 0.0, dyn8.data(), 1);

        std::vector<double> px_rhs_scaled(L);
        std::vector<double> py_rhs_scaled(L);
        std::vector<double> pz_rhs_scaled(L);
        std::vector<double> psi_rhs_scaled(L);

        std::vector<double> vx_rhs_scaled(L);
        std::vector<double> vy_rhs_scaled(L);
        std::vector<double> vz_rhs_scaled(L);
        std::vector<double> w_rhs_scaled(L);

        vdMul(L, vx_vector.data(), r_tau_row_.data(), px_rhs_scaled.data());
        vdMul(L, vy_vector.data(), r_tau_row_.data(), py_rhs_scaled.data());
        vdMul(L, vz_vector.data(), r_tau_row_.data(), pz_rhs_scaled.data());
        vdMul(L, w_vector.data(),  r_tau_row_.data(), psi_rhs_scaled.data());

        vdMul(L, ax_vector.data(), r_tau_row_.data(), vx_rhs_scaled.data());
        vdMul(L, ay_vector.data(), r_tau_row_.data(), vy_rhs_scaled.data());
        vdMul(L, az_vector.data(), r_tau_row_.data(), vz_rhs_scaled.data());
        vdMul(L, aw_vector.data(), r_tau_row_.data(), w_rhs_scaled.data());

        std::vector<double> res_px(L);
        std::vector<double> res_py(L);
        std::vector<double> res_pz(L);
        std::vector<double> res_psi(L);

        std::vector<double> res_vx(L);
        std::vector<double> res_vy(L);
        std::vector<double> res_vz(L);
        std::vector<double> res_w(L);

        vdSub(L, dyn1.data(), px_rhs_scaled.data(),  res_px.data());
        vdSub(L, dyn2.data(), py_rhs_scaled.data(),  res_py.data());
        vdSub(L, dyn3.data(), pz_rhs_scaled.data(),  res_pz.data());
        vdSub(L, dyn4.data(), psi_rhs_scaled.data(), res_psi.data());

        vdSub(L, dyn5.data(), vx_rhs_scaled.data(), res_vx.data());
        vdSub(L, dyn6.data(), vy_rhs_scaled.data(), res_vy.data());
        vdSub(L, dyn7.data(), vz_rhs_scaled.data(), res_vz.data());
        vdSub(L, dyn8.data(), w_rhs_scaled.data(),  res_w.data());

        const double v_max2 = v_max_ * v_max_;
        const double a_max2 = a_max_ * a_max_;

        std::vector<double> vx2(L);
        std::vector<double> vy2(L);
        std::vector<double> vz2(L);

        std::vector<double> ax2(L);
        std::vector<double> ay2(L);
        std::vector<double> az2(L);

        vdSqr(L, vx_vector.data(), vx2.data());
        vdSqr(L, vy_vector.data(), vy2.data());
        vdSqr(L, vz_vector.data(), vz2.data());

        vdSqr(L, ax_vector.data(), ax2.data());
        vdSqr(L, ay_vector.data(), ay2.data());
        vdSqr(L, az_vector.data(), az2.data());

        std::vector<double> v2_sum(L);
        std::vector<double> a2_sum(L);

        vdAdd(L, vx2.data(), vy2.data(), v2_sum.data());
        vdAdd(L, v2_sum.data(), vz2.data(), v2_sum.data());

        vdAdd(L, ax2.data(), ay2.data(), a2_sum.data());
        vdAdd(L, a2_sum.data(), az2.data(), a2_sum.data());

        std::vector<double> v_max2_vec(L, v_max2);
        std::vector<double> a_max2_vec(L, a_max2);

        std::vector<double> c_speed(L);
        std::vector<double> c_accel(L);

        vdSub(L, v2_sum.data(), v_max2_vec.data(), c_speed.data());
        vdSub(L, a2_sum.data(), a_max2_vec.data(), c_accel.data());

        // ------------------------------------------------------------
        // Obstacle avoidance constraints:
        //
        // 1) Cylinder-like obstacle:
        //      (px - cyl_x)^2 + (py - cyl_y)^2 >= cyl_radius^2
        //
        // 2) Sphere obstacle:
        //      (px - sx)^2 + (py - sy)^2 + (pz - sz)^2 >= sphere_radius^2
        //
        // Written for IPOPT as:
        //      c_obs = -dist2_elevated + radius^2 <= 0
        // ------------------------------------------------------------

        auto add_vectors =
            [](const std::vector<double>& a,
               const std::vector<double>& b) -> std::vector<double> {
                const size_t size = std::min(a.size(), b.size());
                std::vector<double> out(size, 0.0);

                for (size_t i = 0; i < size; ++i) {
                    out[i] = a[i] + b[i];
                }

                return out;
            };

        auto binom_ld =
            [](int n, int k) -> long double {
                if (k < 0 || k > n) {
                    return 0.0L;
                }

                if (k == 0 || k == n) {
                    return 1.0L;
                }

                if (k > n - k) {
                    k = n - k;
                }

                long double result = 1.0L;

                for (int i = 1; i <= k; ++i) {
                    result *= static_cast<long double>(n - k + i);
                    result /= static_cast<long double>(i);
                }

                return result;
            };

        auto degree_elevate =
            [&](const std::vector<double>& cp,
                int target_degree) -> std::vector<double> {
                const int current_degree = static_cast<int>(cp.size()) - 1;

                if (target_degree <= current_degree) {
                    return cp;
                }

                const int r = target_degree - current_degree;

                std::vector<double> elevated(target_degree + 1, 0.0);

                for (int j = 0; j <= target_degree; ++j) {
                    const int i_min = std::max(0, j - r);
                    const int i_max = std::min(current_degree, j);

                    long double value = 0.0L;

                    for (int i = i_min; i <= i_max; ++i) {
                        const long double coeff =
                            binom_ld(current_degree, i)
                        * binom_ld(r, j - i)
                        / binom_ld(target_degree, j);

                        value += coeff * static_cast<long double>(cp[i]);
                    }

                    elevated[j] = static_cast<double>(value);
                }

                return elevated;
            };

        std::vector<double> c_obs_all;

        // ------------------------------------------------------------
        // Cylinder obstacle block
        // ------------------------------------------------------------
        {
            std::vector<double> dx(L);
            std::vector<double> dy(L);

            for (int i = 0; i < L; ++i) {
                dx[i] = px_vector[i] - cyl_x_;
                dy[i] = py_vector[i] - cyl_y_;
            }

            std::vector<double> dx2 = BernsteinProduct(dx, dx);
            std::vector<double> dy2 = BernsteinProduct(dy, dy);

            std::vector<double> dist2_cyl =
                add_vectors(dx2, dy2);

            const int deg_current =
                static_cast<int>(dist2_cyl.size()) - 1;

            const int deg_target =
                deg_current + obs_deg_elev_extra_;

            std::vector<double> dist2_cyl_elev =
                degree_elevate(dist2_cyl, deg_target);

            const double cyl_radius2 = cyl_radius_ * cyl_radius_;

            for (double coeff : dist2_cyl_elev) {
                c_obs_all.push_back(-coeff + cyl_radius2);
            }
        }

        // ------------------------------------------------------------
        // Sphere obstacle block
        // ------------------------------------------------------------
        {
            std::vector<double> dx(L);
            std::vector<double> dy(L);
            std::vector<double> dz(L);

            for (int i = 0; i < L; ++i) {
                dx[i] = px_vector[i] - sphere_x_;
                dy[i] = py_vector[i] - sphere_y_;
                dz[i] = pz_vector[i] - sphere_z_;
            }

            std::vector<double> dx2 = BernsteinProduct(dx, dx);
            std::vector<double> dy2 = BernsteinProduct(dy, dy);
            std::vector<double> dz2 = BernsteinProduct(dz, dz);

            std::vector<double> dist2_xy =
                add_vectors(dx2, dy2);

            std::vector<double> dist2_sphere =
                add_vectors(dist2_xy, dz2);

            const int deg_current =
                static_cast<int>(dist2_sphere.size()) - 1;

            const int deg_target =
                deg_current + obs_deg_elev_extra_;

            std::vector<double> dist2_sphere_elev =
                degree_elevate(dist2_sphere, deg_target);

            const double sphere_radius2 = sphere_radius_ * sphere_radius_;

            for (double coeff : dist2_sphere_elev) {
                c_obs_all.push_back(-coeff + sphere_radius2);
            }
        }

        for (Index i = 0; i < L; ++i) {
            g[0 * L + i] = res_px[i];
            g[1 * L + i] = res_py[i];
            g[2 * L + i] = res_pz[i];
            g[3 * L + i] = res_psi[i];

            g[4 * L + i] = res_vx[i];
            g[5 * L + i] = res_vy[i];
            g[6 * L + i] = res_vz[i];
            g[7 * L + i] = res_w[i];

            g[8 * L + i] = c_speed[i];
            g[9 * L + i] = c_accel[i];
        }

        const Index obs_start = 10 * L;

        for (Index i = 0; i < static_cast<Index>(c_obs_all.size()); ++i) {
            g[obs_start + i] = c_obs_all[i];
        }

        return true;
    }

    virtual bool eval_jac_g(
        Index n,
        const Number* x,
        bool new_x,
        Index m,
        Index nele_jac,
        Index* iRow,
        Index* jCol,
        Number* values
    ) {
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

    virtual bool eval_grad_f(
        Index n,
        const Number* x,
        bool new_x,
        Number* grad_f
    ) {
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
        const int L = N_ + 1;

        solution_x_.resize(12 * L);

        for (Index i = 0; i < 12 * L; ++i) {
            solution_x_[i] = x[i];
        }

        final_obj_value_ = obj_value;

        const int K = 1000;

        std::vector<double> t_norm(K);
        std::vector<double> t_real(K);

        for (int i = 0; i < K; ++i) {
            const double s =
                static_cast<double>(i) / static_cast<double>(K - 1);

            t_norm[i] = s;
            t_real[i] = s * tf_;
        }

        std::vector<double> px_vector(
            solution_x_.begin(),
            solution_x_.begin() + L
        );

        std::vector<double> py_vector(
            solution_x_.begin() + 1 * L,
            solution_x_.begin() + 2 * L
        );

        std::vector<double> pz_vector(
            solution_x_.begin() + 2 * L,
            solution_x_.begin() + 3 * L
        );

        std::vector<double> psi_vector(
            solution_x_.begin() + 3 * L,
            solution_x_.begin() + 4 * L
        );

        std::vector<double> vx_vector(
            solution_x_.begin() + 4 * L,
            solution_x_.begin() + 5 * L
        );

        std::vector<double> vy_vector(
            solution_x_.begin() + 5 * L,
            solution_x_.begin() + 6 * L
        );

        std::vector<double> vz_vector(
            solution_x_.begin() + 6 * L,
            solution_x_.begin() + 7 * L
        );

        std::vector<double> w_vector(
            solution_x_.begin() + 7 * L,
            solution_x_.begin() + 8 * L
        );

        std::vector<double> ax_vector(
            solution_x_.begin() + 8 * L,
            solution_x_.begin() + 9 * L
        );

        std::vector<double> ay_vector(
            solution_x_.begin() + 9 * L,
            solution_x_.begin() + 10 * L
        );

        std::vector<double> az_vector(
            solution_x_.begin() + 10 * L,
            solution_x_.begin() + 11 * L
        );

        std::vector<double> aw_vector(
            solution_x_.begin() + 11 * L,
            solution_x_.begin() + 12 * L
        );

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

        auto px_real =
            BernsteinPoly(px_2d, t_real, 0.0, tf_);

        auto py_real =
            BernsteinPoly(py_2d, t_real, 0.0, tf_);

        auto pz_real =
            BernsteinPoly(pz_2d, t_real, 0.0, tf_);

        auto psi_real =
            BernsteinPoly(psi_2d, t_real, 0.0, tf_);

        auto vx_real =
            BernsteinPoly(vx_2d, t_real, 0.0, tf_);

        auto vy_real =
            BernsteinPoly(vy_2d, t_real, 0.0, tf_);

        auto vz_real =
            BernsteinPoly(vz_2d, t_real, 0.0, tf_);

        auto w_real =
            BernsteinPoly(w_2d, t_real, 0.0, tf_);

        auto ax_real =
            BernsteinPoly(ax_2d, t_real, 0.0, tf_);

        auto ay_real =
            BernsteinPoly(ay_2d, t_real, 0.0, tf_);

        auto az_real =
            BernsteinPoly(az_2d, t_real, 0.0, tf_);

        auto aw_real =
            BernsteinPoly(aw_2d, t_real, 0.0, tf_);

        auto flatten =
            [](const std::vector<std::vector<double>>& input) {
                std::vector<double> output;

                for (const auto& row : input) {
                    output.insert(output.end(), row.begin(), row.end());
                }

                return output;
            };

        writeToCSV(t_real, flatten(px_real), "px_real.csv");
        writeToCSV(bebot_.getNodes(), px_vector, "px_controlpoints_real.csv");

        writeToCSV(t_real, flatten(py_real), "py_real.csv");
        writeToCSV(bebot_.getNodes(), py_vector, "py_controlpoints_real.csv");

        writeToCSV(t_real, flatten(pz_real), "pz_real.csv");
        writeToCSV(bebot_.getNodes(), pz_vector, "pz_controlpoints_real.csv");

        writeToCSV(t_real, flatten(psi_real), "psi_real.csv");
        writeToCSV(bebot_.getNodes(), psi_vector, "psi_controlpoints_real.csv");

        writeToCSV(t_real, flatten(vx_real), "vx_real.csv");
        writeToCSV(bebot_.getNodes(), vx_vector, "vx_controlpoints_real.csv");

        writeToCSV(t_real, flatten(vy_real), "vy_real.csv");
        writeToCSV(bebot_.getNodes(), vy_vector, "vy_controlpoints_real.csv");

        writeToCSV(t_real, flatten(vz_real), "vz_real.csv");
        writeToCSV(bebot_.getNodes(), vz_vector, "vz_controlpoints_real.csv");

        writeToCSV(t_real, flatten(w_real), "w_real.csv");
        writeToCSV(bebot_.getNodes(), w_vector, "w_controlpoints_real.csv");

        writeToCSV(t_real, flatten(ax_real), "ax_real.csv");
        writeToCSV(bebot_.getNodes(), ax_vector, "ax_controlpoints_real.csv");

        writeToCSV(t_real, flatten(ay_real), "ay_real.csv");
        writeToCSV(bebot_.getNodes(), ay_vector, "ay_controlpoints_real.csv");

        writeToCSV(t_real, flatten(az_real), "az_real.csv");
        writeToCSV(bebot_.getNodes(), az_vector, "az_controlpoints_real.csv");

        writeToCSV(t_real, flatten(aw_real), "aw_real.csv");
        writeToCSV(bebot_.getNodes(), aw_vector, "aw_controlpoints_real.csv");
    }

    const std::vector<Number>& get_solution_x() const {
        return solution_x_;
    }

    Number get_final_obj_value() const {
        return final_obj_value_;
    }

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

    double cyl_x_;
    double cyl_y_;
    double cyl_radius_;

    double sphere_x_;
    double sphere_y_;
    double sphere_z_;
    double sphere_radius_;

    int obs_deg_elev_extra_;

    Bebot bebot_;

    std::vector<double> tau_vec_;
    std::vector<double> r_tau_;
    std::vector<double> r_tau_row_;

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
        return bernsteinpoly_resultz_;
    }
};

extern "C" {
    PointSetProblem* create_point_set_problem(
        int N,
        double tf,
        double px_max,
        double px_min,
        double py_max,
        double py_min,
        double pz_max,
        double pz_min,
        double psi_max,
        double psi_min,
        double v_max,
        double w_max,
        double a_max,
        double aw_max,
        double px_cur,
        double py_cur,
        double pz_cur,
        double psi_cur,
        double vx_cur,
        double vy_cur,
        double vz_cur,
        double w_cur,
        double pxf,
        double pyf,
        double pzf,
        double psif,

        double cyl_x,
        double cyl_y,
        double cyl_radius,

        double sphere_x,
        double sphere_y,
        double sphere_z,
        double sphere_radius
    ) {
        PointSetProblem* problem = new PointSetProblem(
            N,
            tf,
            px_max,
            px_min,
            py_max,
            py_min,
            pz_max,
            pz_min,
            psi_max,
            psi_min,
            v_max,
            w_max,
            a_max,
            aw_max,
            px_cur,
            py_cur,
            pz_cur,
            psi_cur,
            vx_cur,
            vy_cur,
            vz_cur,
            w_cur,
            pxf,
            pyf,
            pzf,
            psif,

            cyl_x,
            cyl_y,
            cyl_radius,

            sphere_x,
            sphere_y,
            sphere_z,
            sphere_radius
        );

        /*
            IMPORTANT LIFETIME FIX:

            IPOPT uses SmartPtr reference counting for TNLP objects.
            The external C caller receives a raw pointer, so we keep one
            explicit reference for that external owner.

            solve_point_set_problem() will create a temporary SmartPtr<TNLP>,
            which adds/releases its own reference. This external reference
            keeps the object alive after OptimizeTNLP() returns.

            destroy_point_set_problem() releases this external reference.
        */
        problem->AddRef(nullptr);

        return problem;
    }

    void solve_point_set_problem(PointSetProblem* problem) {
        if (!problem) {
            std::cerr << "solve_point_set_problem() received nullptr"
                      << std::endl;
            return;
        }

        SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

        app->Options()->SetStringValue("linear_solver", "ma57");
        app->Options()->SetStringValue("mu_strategy", "adaptive");

        app->Options()->SetStringValue(
            "gradient_approximation",
            "finite-difference-values"
        );

        app->Options()->SetStringValue(
            "jacobian_approximation",
            "finite-difference-values"
        );

        app->Options()->SetStringValue(
            "hessian_approximation",
            "limited-memory"
        );

        app->Options()->SetIntegerValue("max_iter", 400);

        app->Options()->SetNumericValue("tol", 1e-3);
        app->Options()->SetNumericValue("constr_viol_tol", 1e-3);
        app->Options()->SetNumericValue("obj_scaling_factor", 1e-3);

        app->Options()->SetIntegerValue("print_level", 3);

        app->RethrowNonIpoptException(true);

        ApplicationReturnStatus status = app->Initialize();

        if (status != Solve_Succeeded) {
            std::cerr << "IPOPT initialization failed!" << std::endl;
            return;
        }

        SmartPtr<TNLP> tnlp = problem;
        status = app->OptimizeTNLP(tnlp);

        if (status == Solve_Succeeded || status == Solved_To_Acceptable_Level) {
            std::cout << "Optimization succeeded!" << std::endl;

            const auto& solution_x = problem->get_solution_x();

            std::cout << "Optimal Solution (x): ";

            for (Index i = 0; i < static_cast<Index>(solution_x.size()); i++) {
                std::cout << solution_x[i] << " ";
            }

            std::cout << std::endl;
        } else {
            std::cerr << "Optimization failed with status "
                      << status << std::endl;
        }
    }

    void get_solution(PointSetProblem* problem, double* solution, int n) {
        if (!problem) {
            std::cerr << "get_solution() received nullptr problem"
                      << std::endl;
            return;
        }

        if (!solution) {
            std::cerr << "get_solution() received nullptr solution buffer"
                      << std::endl;
            return;
        }

        if (n <= 0) {
            std::cerr << "get_solution() received non-positive n = "
                      << n << std::endl;
            return;
        }

        const std::vector<double>& sol = problem->get_solution_x();

        std::cout << "[DEBUG] get_solution called. Vector size: "
                  << sol.size() << ", requested n: " << n << std::endl;

        const int copy_n =
            std::min(n, static_cast<int>(sol.size()));

        std::copy(sol.begin(), sol.begin() + copy_n, solution);

        if (copy_n < n) {
            std::fill(solution + copy_n, solution + n, 0.0);

            std::cerr << "[WARN] get_solution requested " << n
                      << " values, but solution has only " << sol.size()
                      << ". Remaining output entries filled with 0."
                      << std::endl;
        }
    }

    double get_final_objective_value(PointSetProblem* problem) {
        if (!problem) {
            std::cerr << "get_final_objective_value() received nullptr"
                      << std::endl;
            return std::numeric_limits<double>::quiet_NaN();
        }

        return problem->get_final_obj_value();
    }

    void destroy_point_set_problem(PointSetProblem* problem) {
        if (problem) {
            problem->ReleaseRef(nullptr);
        }
    }
}