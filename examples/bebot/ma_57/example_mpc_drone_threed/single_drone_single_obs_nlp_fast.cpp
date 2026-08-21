#include "../../../../Ipopt_ma57_solver/src/Interfaces/IpIpoptApplication.hpp"
#include "../../../../Ipopt_ma57_solver/src/Interfaces/IpTNLP.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include "../../../../include/bebot.h"

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

        double sphere_x,
        double sphere_y,
        double sphere_z,
        double sphere_radius
    )
        : N_(N),
          L_(N + 1),
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

          sphere_x_(sphere_x),
          sphere_y_(sphere_y),
          sphere_z_(sphere_z),
          sphere_radius_(sphere_radius),

          obs_deg_elev_extra_(5),
          product_degree_(2 * N),
          sphere_degree_(2 * N + obs_deg_elev_extra_),
          sphere_constraint_count_(sphere_degree_ + 1),

          bebot_(N, tf_),

          final_obj_value_(
              std::numeric_limits<double>::quiet_NaN()
          ) {

        bebot_.calculate();

        // --------------------------------------------------------
        // Time-transformation coefficients.
        // Same mathematical definition as the original code.
        // --------------------------------------------------------
        r_tau_.resize(L_);

        if (L_ == 1) {
            const double tau = -1.0;
            const double denom = 1.0 - tau;
            r_tau_[0] = 2.0 / (denom * denom);
        } else {
            const double tau0 = -1.0;
            const double tau1 = 0.9;
            const double step =
                (tau1 - tau0)
                / static_cast<double>(N_);

            for (int i = 0; i < L_; ++i) {
                const double tau =
                    tau0
                    + step * static_cast<double>(i);

                const double denom =
                    1.0 - tau;

                r_tau_[i] =
                    2.0 / (denom * denom);
            }
        }

        // --------------------------------------------------------
        // Precompute Bernstein-product weights.
        //
        // For two degree-N Bernstein polynomials:
        //
        // q_k =
        //   sum_i [
        //      C(N,i) C(N,k-i) / C(2N,k)
        //   ] a_i b_{k-i}.
        //
        // Stored as:
        // product_weights_[k * L_ + i].
        // --------------------------------------------------------
        product_weights_.assign(
            (product_degree_ + 1) * L_,
            0.0
        );

        for (int k = 0; k <= product_degree_; ++k) {
            const int i_min =
                std::max(0, k - N_);

            const int i_max =
                std::min(N_, k);

            const long double denom =
                binom_ld(product_degree_, k);

            for (int i = i_min; i <= i_max; ++i) {
                const int other =
                    k - i;

                const long double weight =
                    binom_ld(N_, i)
                    * binom_ld(N_, other)
                    / denom;

                product_weights_[
                    k * L_ + i
                ] =
                    static_cast<double>(weight);
            }
        }

        // --------------------------------------------------------
        // Precompute degree-elevation matrix.
        //
        // Degree product_degree_ -> sphere_degree_.
        //
        // e_j = sum_k E(j,k) q_k
        //
        // Stored as:
        // elevation_weights_[j * (product_degree_ + 1) + k].
        // --------------------------------------------------------
        const int elevation_amount =
            sphere_degree_ - product_degree_;

        elevation_weights_.assign(
            sphere_constraint_count_
                * (product_degree_ + 1),
            0.0
        );

        for (
            int j = 0;
            j <= sphere_degree_;
            ++j
        ) {
            const long double denom =
                binom_ld(sphere_degree_, j);

            const int k_min =
                std::max(
                    0,
                    j - elevation_amount
                );

            const int k_max =
                std::min(
                    product_degree_,
                    j
                );

            for (
                int k = k_min;
                k <= k_max;
                ++k
            ) {
                const long double weight =
                    binom_ld(product_degree_, k)
                    * binom_ld(
                        elevation_amount,
                        j - k
                    )
                    / denom;

                elevation_weights_[
                    j * (product_degree_ + 1)
                    + k
                ] =
                    static_cast<double>(weight);
            }
        }

        // Scratch storage reused by eval_g().
        sphere_product_coeffs_.resize(
            product_degree_ + 1,
            0.0
        );
    }

    // ========================================================
    // NLP dimensions / sparsity
    // ========================================================
    bool get_nlp_info(
        Index& n,
        Index& m,
        Index& nnz_jac_g,
        Index& nnz_h_lag,
        IndexStyleEnum& index_style
    ) override {
        n =
            12 * L_;

        m =
            10 * L_
            + sphere_constraint_count_;

        // ----------------------------------------------------
        // Exact sparse Jacobian structure.
        //
        // Dynamics:
        // 8*L rows.
        // Each row depends on:
        //   L coefficients of differentiated state
        //   1 matching velocity/control coefficient
        //
        // Speed:
        // 3 entries per row.
        //
        // Acceleration:
        // 3 entries per row.
        //
        // Sphere:
        // each elevated sphere coefficient depends on
        // all px, py, pz Bernstein coefficients = 3L.
        // ----------------------------------------------------
        const Index dynamics_nnz =
            8 * L_ * (L_ + 1);

        const Index speed_nnz =
            3 * L_;

        const Index acceleration_nnz =
            3 * L_;

        const Index sphere_nnz =
            sphere_constraint_count_
            * 3 * L_;

        nnz_jac_g =
            dynamics_nnz
            + speed_nnz
            + acceleration_nnz
            + sphere_nnz;

        // Limited-memory Hessian.
        nnz_h_lag = 0;

        index_style =
            TNLP::C_STYLE;

        return true;
    }

    // ========================================================
    // Variable and constraint bounds
    // ========================================================
    bool get_bounds_info(
        Index n,
        Number* x_l,
        Number* x_u,
        Index m,
        Number* g_l,
        Number* g_u
    ) override {
        const double inf =
            std::numeric_limits<double>::infinity();

        for (Index i = 0; i < n; ++i) {
            x_l[i] = -inf;
            x_u[i] =  inf;
        }

        // ----------------------------------------------------
        // px
        // ----------------------------------------------------
        x_l[0 * L_] = px_cur_;
        x_u[0 * L_] = px_cur_;

        for (int i = 1; i < L_; ++i) {
            x_l[0 * L_ + i] = px_min_;
            x_u[0 * L_ + i] = px_max_;
        }

        // ----------------------------------------------------
        // py
        // ----------------------------------------------------
        x_l[1 * L_] = py_cur_;
        x_u[1 * L_] = py_cur_;

        for (int i = 1; i < L_; ++i) {
            x_l[1 * L_ + i] = py_min_;
            x_u[1 * L_ + i] = py_max_;
        }

        // ----------------------------------------------------
        // pz
        // ----------------------------------------------------
        x_l[2 * L_] = pz_cur_;
        x_u[2 * L_] = pz_cur_;

        for (int i = 1; i < L_; ++i) {
            x_l[2 * L_ + i] = pz_min_;
            x_u[2 * L_ + i] = pz_max_;
        }

        // ----------------------------------------------------
        // psi
        // ----------------------------------------------------
        x_l[3 * L_] = psi_cur_;
        x_u[3 * L_] = psi_cur_;

        for (int i = 1; i < L_; ++i) {
            x_l[3 * L_ + i] = psi_min_;
            x_u[3 * L_ + i] = psi_max_;
        }

        // ----------------------------------------------------
        // vx, vy, vz
        // ----------------------------------------------------
        const double current_v[3] = {
            vx_cur_,
            vy_cur_,
            vz_cur_
        };

        for (int block = 4; block <= 6; ++block) {
            x_l[block * L_] =
                current_v[block - 4];

            x_u[block * L_] =
                current_v[block - 4];

            for (int i = 1; i < L_; ++i) {
                x_l[block * L_ + i] =
                    -v_max_;

                x_u[block * L_ + i] =
                     v_max_;
            }
        }

        // ----------------------------------------------------
        // yaw rate w
        // ----------------------------------------------------
        x_l[7 * L_] = w_cur_;
        x_u[7 * L_] = w_cur_;

        for (int i = 1; i < L_; ++i) {
            x_l[7 * L_ + i] =
                -w_max_;

            x_u[7 * L_ + i] =
                 w_max_;
        }

        // ----------------------------------------------------
        // ax, ay, az
        // ----------------------------------------------------
        for (int block = 8; block <= 10; ++block) {
            for (int i = 0; i < L_; ++i) {
                x_l[block * L_ + i] =
                    -a_max_;

                x_u[block * L_ + i] =
                     a_max_;
            }
        }

        // ----------------------------------------------------
        // aw
        // ----------------------------------------------------
        for (int i = 0; i < L_; ++i) {
            x_l[11 * L_ + i] =
                -aw_max_;

            x_u[11 * L_ + i] =
                 aw_max_;
        }

        // ----------------------------------------------------
        // Dynamics equality constraints.
        // ----------------------------------------------------
        for (int i = 0; i < 8 * L_; ++i) {
            g_l[i] = 0.0;
            g_u[i] = 0.0;
        }

        // ----------------------------------------------------
        // Speed, acceleration, sphere:
        // all written c(x) <= 0.
        // ----------------------------------------------------
        for (
            Index i = 8 * L_;
            i < m;
            ++i
        ) {
            g_l[i] = -inf;
            g_u[i] = 0.0;
        }

        return true;
    }

    // ========================================================
    // Starting point
    // ========================================================
    bool get_starting_point(
        Index n,
        bool init_x,
        Number* x,
        bool init_z,
        Number* z_L,
        Number* z_U,
        Index m,
        bool init_lambda,
        Number* lambda
    ) override {
        // Position / yaw:
        // straight interpolation current -> target.
        const double q0[4] = {
            px_cur_,
            py_cur_,
            pz_cur_,
            psi_cur_
        };

        const double qf[4] = {
            pxf_,
            pyf_,
            pzf_,
            psif_
        };

        for (int block = 0; block < 4; ++block) {
            for (int i = 0; i < L_; ++i) {
                const double alpha =
                    (N_ == 0)
                        ? 0.0
                        : static_cast<double>(i)
                          / static_cast<double>(N_);

                x[block * L_ + i] =
                    q0[block]
                    + alpha
                    * (qf[block] - q0[block]);
            }
        }

        // Velocity / yaw rate:
        // constant current values.
        const double v0[4] = {
            vx_cur_,
            vy_cur_,
            vz_cur_,
            w_cur_
        };

        for (int block = 4; block < 8; ++block) {
            const double value =
                v0[block - 4];

            for (int i = 0; i < L_; ++i) {
                x[block * L_ + i] =
                    value;
            }
        }

        // Acceleration / yaw acceleration:
        // zero.
        for (int block = 8; block < 12; ++block) {
            for (int i = 0; i < L_; ++i) {
                x[block * L_ + i] =
                    0.0;
            }
        }

        return true;
    }

    // ========================================================
    // Objective
    // ========================================================
    bool eval_f(
        Index n,
        const Number* x,
        bool new_x,
        Number& obj_value
    ) override {
        constexpr double w_p =
            0.01;

        constexpr double w_psi =
            1.0;

        constexpr double w_a =
            0.01;

        constexpr double w_aw =
            0.001;

        obj_value = 0.0;

        // Position objective.
        for (int i = 0; i < L_; ++i) {
            const double dx =
                x[0 * L_ + i] - pxf_;

            const double dy =
                x[1 * L_ + i] - pyf_;

            const double dz =
                x[2 * L_ + i] - pzf_;

            const double dpsi =
                x[3 * L_ + i] - psif_;

            obj_value +=
                r_tau_[i]
                * (
                    w_p
                    * (
                        dx * dx
                        + dy * dy
                        + dz * dz
                    )
                    + w_psi
                    * dpsi * dpsi
                );
        }

        // Translational acceleration objective.
        for (int i = 0; i < L_; ++i) {
            const double ax =
                x[8 * L_ + i];

            const double ay =
                x[9 * L_ + i];

            const double az =
                x[10 * L_ + i];

            const double aw =
                x[11 * L_ + i];

            obj_value +=
                w_a
                * (
                    ax * ax
                    + ay * ay
                    + az * az
                )
                + w_aw
                * aw * aw;
        }

        return true;
    }

    // ========================================================
    // Exact objective gradient
    // ========================================================
    bool eval_grad_f(
        Index n,
        const Number* x,
        bool new_x,
        Number* grad_f
    ) override {
        constexpr double w_p =
            0.01;

        constexpr double w_psi =
            1.0;

        constexpr double w_a =
            0.01;

        constexpr double w_aw =
            0.001;

        std::fill(
            grad_f,
            grad_f + n,
            0.0
        );

        for (int i = 0; i < L_; ++i) {
            grad_f[0 * L_ + i] =
                2.0
                * w_p
                * r_tau_[i]
                * (
                    x[0 * L_ + i]
                    - pxf_
                );

            grad_f[1 * L_ + i] =
                2.0
                * w_p
                * r_tau_[i]
                * (
                    x[1 * L_ + i]
                    - pyf_
                );

            grad_f[2 * L_ + i] =
                2.0
                * w_p
                * r_tau_[i]
                * (
                    x[2 * L_ + i]
                    - pzf_
                );

            grad_f[3 * L_ + i] =
                2.0
                * w_psi
                * r_tau_[i]
                * (
                    x[3 * L_ + i]
                    - psif_
                );

            grad_f[8 * L_ + i] =
                2.0
                * w_a
                * x[8 * L_ + i];

            grad_f[9 * L_ + i] =
                2.0
                * w_a
                * x[9 * L_ + i];

            grad_f[10 * L_ + i] =
                2.0
                * w_a
                * x[10 * L_ + i];

            grad_f[11 * L_ + i] =
                2.0
                * w_aw
                * x[11 * L_ + i];
        }

        return true;
    }

    // ========================================================
    // Constraints
    // ========================================================
    bool eval_g(
        Index n,
        const Number* x,
        bool new_x,
        Index m,
        Number* g
    ) override {
        const auto& Dm =
            bebot_.getDifferentiationMatrix();

        // ----------------------------------------------------
        // 8 dynamics blocks.
        //
        // Existing code used:
        //
        // cblas_dgemv(
        //     CblasColMajor,
        //     CblasTrans,
        //     ...
        // )
        //
        // Therefore:
        //
        // (Dm^T q)_i =
        //     sum_j Dm[j + i*L] q_j.
        //
        // This direct loop preserves that convention.
        // ----------------------------------------------------
        for (int block = 0; block < 8; ++block) {
            const double* q =
                x + block * L_;

            const double* rhs =
                x + (block + 4) * L_;

            for (int i = 0; i < L_; ++i) {
                double derivative =
                    0.0;

                const int col_offset =
                    i * L_;

                for (int j = 0; j < L_; ++j) {
                    derivative +=
                        Dm[col_offset + j]
                        * q[j];
                }

                g[block * L_ + i] =
                    derivative
                    - r_tau_[i]
                    * rhs[i];
            }
        }

        // ----------------------------------------------------
        // Translational speed:
        // vx^2 + vy^2 + vz^2 - v_max^2 <= 0.
        // ----------------------------------------------------
        const double v_max2 =
            v_max_ * v_max_;

        for (int i = 0; i < L_; ++i) {
            const double vx =
                x[4 * L_ + i];

            const double vy =
                x[5 * L_ + i];

            const double vz =
                x[6 * L_ + i];

            g[8 * L_ + i] =
                vx * vx
                + vy * vy
                + vz * vz
                - v_max2;
        }

        // ----------------------------------------------------
        // Translational acceleration:
        // ax^2 + ay^2 + az^2 - a_max^2 <= 0.
        // ----------------------------------------------------
        const double a_max2 =
            a_max_ * a_max_;

        for (int i = 0; i < L_; ++i) {
            const double ax =
                x[8 * L_ + i];

            const double ay =
                x[9 * L_ + i];

            const double az =
                x[10 * L_ + i];

            g[9 * L_ + i] =
                ax * ax
                + ay * ay
                + az * az
                - a_max2;
        }

        // ----------------------------------------------------
        // One 3-D spherical obstacle.
        //
        // q_k = Bernstein coefficients of
        //
        //   (px-sx)^2
        // + (py-sy)^2
        // + (pz-sz)^2.
        //
        // The Bernstein product weights and degree-elevation
        // matrix were precomputed once in the constructor.
        // ----------------------------------------------------
        for (
            int k = 0;
            k <= product_degree_;
            ++k
        ) {
            double qk =
                0.0;

            const int i_min =
                std::max(0, k - N_);

            const int i_max =
                std::min(N_, k);

            for (int i = i_min; i <= i_max; ++i) {
                const int other =
                    k - i;

                const double weight =
                    product_weights_[
                        k * L_ + i
                    ];

                const double dx_i =
                    x[0 * L_ + i]
                    - sphere_x_;

                const double dx_other =
                    x[0 * L_ + other]
                    - sphere_x_;

                const double dy_i =
                    x[1 * L_ + i]
                    - sphere_y_;

                const double dy_other =
                    x[1 * L_ + other]
                    - sphere_y_;

                const double dz_i =
                    x[2 * L_ + i]
                    - sphere_z_;

                const double dz_other =
                    x[2 * L_ + other]
                    - sphere_z_;

                qk +=
                    weight
                    * (
                        dx_i * dx_other
                        + dy_i * dy_other
                        + dz_i * dz_other
                    );
            }

            sphere_product_coeffs_[k] =
                qk;
        }

        const double radius2 =
            sphere_radius_
            * sphere_radius_;

        const Index sphere_start =
            10 * L_;

        for (
            int j = 0;
            j < sphere_constraint_count_;
            ++j
        ) {
            double elevated =
                0.0;

            const int elevation_offset =
                j * (product_degree_ + 1);

            for (
                int k = 0;
                k <= product_degree_;
                ++k
            ) {
                elevated +=
                    elevation_weights_[
                        elevation_offset + k
                    ]
                    * sphere_product_coeffs_[k];
            }

            g[sphere_start + j] =
                radius2
                - elevated;
        }

        return true;
    }

    // ========================================================
    // Exact sparse constraint Jacobian
    // ========================================================
    bool eval_jac_g(
        Index n,
        const Number* x,
        bool new_x,
        Index m,
        Index nele_jac,
        Index* iRow,
        Index* jCol,
        Number* values
    ) override {
        Index idx = 0;

        // ----------------------------------------------------
        // Structure request.
        // ----------------------------------------------------
        if (values == nullptr) {
            // Dynamics.
            for (int block = 0; block < 8; ++block) {
                for (int i = 0; i < L_; ++i) {
                    const Index row =
                        block * L_ + i;

                    // Dm^T contribution:
                    // all coefficients of q block.
                    for (int j = 0; j < L_; ++j) {
                        iRow[idx] = row;
                        jCol[idx] =
                            block * L_ + j;
                        ++idx;
                    }

                    // Matching rhs coefficient.
                    iRow[idx] = row;
                    jCol[idx] =
                        (block + 4) * L_ + i;
                    ++idx;
                }
            }

            // Speed.
            for (int i = 0; i < L_; ++i) {
                const Index row =
                    8 * L_ + i;

                for (int block = 4; block <= 6; ++block) {
                    iRow[idx] = row;
                    jCol[idx] =
                        block * L_ + i;
                    ++idx;
                }
            }

            // Acceleration.
            for (int i = 0; i < L_; ++i) {
                const Index row =
                    9 * L_ + i;

                for (int block = 8; block <= 10; ++block) {
                    iRow[idx] = row;
                    jCol[idx] =
                        block * L_ + i;
                    ++idx;
                }
            }

            // Sphere.
            for (
                int j = 0;
                j < sphere_constraint_count_;
                ++j
            ) {
                const Index row =
                    10 * L_ + j;

                for (int block = 0; block < 3; ++block) {
                    for (int i = 0; i < L_; ++i) {
                        iRow[idx] = row;
                        jCol[idx] =
                            block * L_ + i;
                        ++idx;
                    }
                }
            }

            return idx == nele_jac;
        }

        // ----------------------------------------------------
        // Numerical values.
        // ----------------------------------------------------
        const auto& Dm =
            bebot_.getDifferentiationMatrix();

        // Dynamics.
        for (int block = 0; block < 8; ++block) {
            for (int i = 0; i < L_; ++i) {
                const int col_offset =
                    i * L_;

                for (int j = 0; j < L_; ++j) {
                    values[idx++] =
                        Dm[col_offset + j];
                }

                values[idx++] =
                    -r_tau_[i];
            }
        }

        // Speed.
        for (int i = 0; i < L_; ++i) {
            values[idx++] =
                2.0 * x[4 * L_ + i];

            values[idx++] =
                2.0 * x[5 * L_ + i];

            values[idx++] =
                2.0 * x[6 * L_ + i];
        }

        // Acceleration.
        for (int i = 0; i < L_; ++i) {
            values[idx++] =
                2.0 * x[8 * L_ + i];

            values[idx++] =
                2.0 * x[9 * L_ + i];

            values[idx++] =
                2.0 * x[10 * L_ + i];
        }

        // ----------------------------------------------------
        // Sphere Jacobian.
        //
        // For one coordinate d_i = p_i - center:
        //
        // q_k =
        //   sum_i w(k,i) d_i d_{k-i}
        //
        // Because w(k,i) is symmetric:
        //
        // dq_k / dd_l =
        //   2 w(k,l) d_{k-l}
        //
        // whenever k-l is a valid degree-N index.
        //
        // c_j = r^2 - sum_k E(j,k) q_k
        //
        // dc_j / dp_l =
        //   -sum_k E(j,k) dq_k/dd_l.
        // ----------------------------------------------------
        const double centers[3] = {
            sphere_x_,
            sphere_y_,
            sphere_z_
        };

        for (
            int j = 0;
            j < sphere_constraint_count_;
            ++j
        ) {
            const int elevation_offset =
                j * (product_degree_ + 1);

            for (int block = 0; block < 3; ++block) {
                const double center =
                    centers[block];

                for (int l = 0; l < L_; ++l) {
                    double derivative_elevated =
                        0.0;

                    const int k_min =
                        l;

                    const int k_max =
                        std::min(
                            product_degree_,
                            l + N_
                        );

                    for (
                        int k = k_min;
                        k <= k_max;
                        ++k
                    ) {
                        const int other =
                            k - l;

                        const double product_weight =
                            product_weights_[
                                k * L_ + l
                            ];

                        const double elevation_weight =
                            elevation_weights_[
                                elevation_offset + k
                            ];

                        const double d_other =
                            x[block * L_ + other]
                            - center;

                        derivative_elevated +=
                            elevation_weight
                            * 2.0
                            * product_weight
                            * d_other;
                    }

                    values[idx++] =
                        -derivative_elevated;
                }
            }
        }

        return idx == nele_jac;
    }

    // ========================================================
    // Final solution
    // ========================================================
    void finalize_solution(
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
    ) override {
        solution_x_.assign(
            x,
            x + 12 * L_
        );

        final_obj_value_ =
            obj_value;
    }

    const std::vector<Number>& get_solution_x() const {
        return solution_x_;
    }

    Number get_final_obj_value() const {
        return final_obj_value_;
    }

private:
    static long double binom_ld(
        int n,
        int k
    ) {
        if (k < 0 || k > n) {
            return 0.0L;
        }

        if (k == 0 || k == n) {
            return 1.0L;
        }

        if (k > n - k) {
            k = n - k;
        }

        long double result =
            1.0L;

        for (int i = 1; i <= k; ++i) {
            result *=
                static_cast<long double>(
                    n - k + i
                );

            result /=
                static_cast<long double>(i);
        }

        return result;
    }

    int N_;
    int L_;

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

    double sphere_x_;
    double sphere_y_;
    double sphere_z_;
    double sphere_radius_;

    int obs_deg_elev_extra_;
    int product_degree_;
    int sphere_degree_;
    int sphere_constraint_count_;

    Bebot bebot_;

    std::vector<double> r_tau_;

    // Precomputed math for sphere constraint.
    std::vector<double> product_weights_;
    std::vector<double> elevation_weights_;

    // Reused scratch storage.
    std::vector<double> sphere_product_coeffs_;

    std::vector<Number> solution_x_;
    Number final_obj_value_;
};

// ============================================================
// C API
// ============================================================
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

    double sphere_x,
    double sphere_y,
    double sphere_z,
    double sphere_radius
) {
    PointSetProblem* problem =
        new PointSetProblem(
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

            sphere_x,
            sphere_y,
            sphere_z,
            sphere_radius
        );

    // Keep one explicit reference for the external raw-pointer owner.
    problem->AddRef(nullptr);

    return problem;
}

int solve_point_set_problem(
    PointSetProblem* problem
) {
    if (!problem) {
        std::cerr
            << "solve_point_set_problem(): nullptr problem\n";

        return 0;
    }

    SmartPtr<IpoptApplication> app =
        IpoptApplicationFactory();

    app->Options()->SetStringValue(
        "linear_solver",
        "ma57"
    );

    app->Options()->SetStringValue(
        "mu_strategy",
        "adaptive"
    );

    // Exact first derivatives are now supplied by this TNLP.
    app->Options()->SetStringValue(
        "gradient_approximation",
        "exact"
    );

    app->Options()->SetStringValue(
        "jacobian_approximation",
        "exact"
    );

    // Keep limited-memory Hessian for now.
    app->Options()->SetStringValue(
        "hessian_approximation",
        "limited-memory"
    );

    app->Options()->SetIntegerValue(
        "max_iter",
        400
    );

    app->Options()->SetNumericValue(
        "tol",
        1e-3
    );

    app->Options()->SetNumericValue(
        "constr_viol_tol",
        1e-3
    );

    app->Options()->SetNumericValue(
        "obj_scaling_factor",
        1e-3
    );

    app->Options()->SetIntegerValue(
        "print_level",
        3
    );

    app->RethrowNonIpoptException(true);

    ApplicationReturnStatus status =
        app->Initialize();

    if (status != Solve_Succeeded) {
        std::cerr
            << "IPOPT initialization failed. Status = "
            << static_cast<int>(status)
            << "\n";

        return 0;
    }

    SmartPtr<TNLP> tnlp =
        problem;

    status =
        app->OptimizeTNLP(tnlp);

    const bool ok =
        status == Solve_Succeeded
        || status == Solved_To_Acceptable_Level;

    if (ok) {
        std::cout
            << "Optimization succeeded. IPOPT status = "
            << static_cast<int>(status)
            << "\n";
    } else {
        std::cerr
            << "Optimization failed. IPOPT status = "
            << static_cast<int>(status)
            << "\n";
    }

    return ok ? 1 : 0;
}

int get_solution_size(
    PointSetProblem* problem
) {
    if (!problem) {
        return 0;
    }

    return static_cast<int>(
        problem->get_solution_x().size()
    );
}

int get_solution(
    PointSetProblem* problem,
    double* solution,
    int n
) {
    if (
        !problem
        || !solution
        || n <= 0
    ) {
        return 0;
    }

    const std::vector<double>& sol =
        problem->get_solution_x();

    const int copy_n =
        std::min(
            n,
            static_cast<int>(sol.size())
        );

    std::copy(
        sol.begin(),
        sol.begin() + copy_n,
        solution
    );

    return copy_n;
}

double get_final_objective_value(
    PointSetProblem* problem
) {
    if (!problem) {
        return
            std::numeric_limits<double>::quiet_NaN();
    }

    return
        problem->get_final_obj_value();
}

void destroy_point_set_problem(
    PointSetProblem* problem
) {
    if (problem) {
        problem->ReleaseRef(nullptr);
    }
}

} // extern "C"
