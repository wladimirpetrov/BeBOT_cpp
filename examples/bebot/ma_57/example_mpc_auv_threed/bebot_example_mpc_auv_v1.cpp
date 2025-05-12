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
    PointSetProblem(int N, double tf, double delta_v_max, double delta_v_min, double delta_s_max, double delta_s_min, double delta_m_max, double delta_m_min, double zmax, double zmin, double wmax, double wmin, double thetamax, double thetamin, double qmax, double qmin, double z0, double w0, double theta0, double q0, double delta_v0, double delta_s0, double delta_m0, double zf, double thetaf, const std::vector<double>& A_flat, const std::vector<double>& B_flat)
        : N_(N), tf_(tf), delta_v_max_(delta_v_max), delta_v_min_(delta_v_min), delta_s_max_(delta_s_max), delta_s_min_(delta_s_min), delta_m_max_(delta_m_max), delta_m_min_(delta_m_min), zmax_(zmax), zmin_(zmin), wmax_(wmax), wmin_(wmin), thetamax_(thetamax), thetamin_(thetamin), qmax_(qmax), qmin_(qmin), z0_(z0), w0_(w0), theta0_(theta0), q0_(q0), delta_v0_(delta_v0), delta_s0_(delta_s0), delta_m0_(delta_m0), zf_(zf), thetaf_(thetaf), bebot_(N, tf_), A_(A_flat), B_(B_flat) {
        bebot_.calculate();
    } 

    // write data to CSV
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
        n = 7 * (N_ + 1); 
        m = 4 * (N_ + 1);
        nnz_jac_g = n*m;  
        nnz_h_lag = 0; 
        index_style = TNLP::C_STYLE;
        return true;
    }

    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u) {
        // x vector (corresponding to z,w,theta,q)
        for (int i = 0; i < 7 * (N_+ 1); i++) {
            x_l[i] = -std::numeric_limits<double>::infinity();
            x_u[i] = std::numeric_limits<double>::infinity();
        }

        // lower and upper bounds for z values excluding the 1st and the last ones
        for (int i = 1; i < N_; ++i) {
            x_l[i] = zmin_;
            x_u[i] = zmax_;
        }
        // lower and upper bounds for the 1st and the last z values
        x_l[0] = z0_;
        x_u[0] = z0_;
        x_l[N_] = zf_;
        x_u[N_] = zf_;

        // lower and upper bounds for w values excluding the 1st one
        for (int i = N_ + 2; i < 2 * (N_ + 1); ++i) {
            x_l[i] = wmin_;
            x_u[i] = wmax_;
        }

        // lower and upper bounds for the 1st w value
        x_l[N_ + 1] = w0_; 
        x_u[N_ + 1] = w0_;

        // lower and upper bounds for theta values excluding the 1st and the last ones
        for (int i = 2 * (N_ + 1) + 1; i < 3 * (N_ + 1) - 1; ++i) {
            x_l[i] = thetamin_;
            x_u[i] = thetamax_;
        }

        // lower and upper bounds for the 1st and the last theta values
        x_l[2 * (N_ + 1)] = theta0_; 
        x_u[2 * (N_ + 1)] = theta0_; 
        x_l[3 * (N_ + 1) - 1] = thetaf_; 
        x_u[3 * (N_ + 1) - 1] = thetaf_;

        // lower and upper bounds for q values excluding the 1st one
        for (int i = 3 * (N_ + 1) + 1; i < 4 * (N_ + 1); ++i) {
            x_l[i] = qmin_;
            x_u[i] = qmax_;
        }

        // lower and upper bounds for the 1st q value
        //x_l[3 * (N_ + 1)] = q0_; 
        //x_u[3 * (N_ + 1)] = q0_;

        // lower and upper bounds for delta_v values excluding the 1st one
        for (int i = 4 * (N_ + 1) + 1; i < 5 * (N_ + 1); ++i) {
            x_l[i] = delta_v_min_;
            x_u[i] = delta_v_max_;
        }
        x_l[4 * (N_ + 1)] = delta_v0_; 
        x_u[4 * (N_ + 1)] = delta_v0_; 

        // lower and upper bounds for delta_s values excluding the 1st one
        for (int i = 5 * (N_ + 1) + 1; i < 6 * (N_ + 1); ++i) {
            x_l[i] = delta_s_min_;
            x_u[i] = delta_s_max_;
        }
        x_l[5 * (N_ + 1)] = delta_s0_; 
        x_u[5 * (N_ + 1)] = delta_s0_; 

        // lower and upper bounds for delta_m values excluding the 1st one
        for (int i = 6 * (N_ + 1) + 1; i < 7 * (N_ + 1); ++i) {
            x_l[i] = delta_m_min_;
            x_u[i] = delta_m_max_;
        }
        x_l[6 * (N_ + 1)] = delta_m0_; 
        x_u[6 * (N_ + 1)] = delta_m0_; 

        // Print constraint bounds if needed
            //for (Index i = 0; i < 7 * (N_ + 1); ++i) {
            //    std::cout << "x_l[" << i << "] = " << x_l[i] << ", x_u[" << i << "] = " << x_u[i] << std::endl;
            //}
        
        // g vector (corresponding to z_dot,w_dot,theta_dot,q_dot,delta_v,delta_s,delta_m)
        
        // first 4 (N_ + 1) vectors z_dot,w_dot,theta_dot,q_dot: -inf to +inf 
        for (int i = 0; i < 4 * (N_ + 1); ++i) {
            g_l[i] = 0;
            g_u[i] = 0;
        }

        // Print constraint bounds if needed
        //for (Index i = 0; i < 4 * (N_ + 1); ++i) {
        //    std::cout << "g_l[" << i << "] = " << g_l[i] << ", g_u[" << i << "] = " << g_u[i] << std::endl;
        //}

        return true;
    }

    // initialization of the starting point
    virtual bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda) {
        for (int i = 0; i < 7 * (N_ + 1); i++) {
            x[i] = 1;
            //std::cout << "x[" << i << "] = " << x[i] << std::endl;
        }
        return true;
    }

    // Objective function
    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value) {
        double w1 = 20.0;
        double w2 = 20.0;
        obj_value = 0.0;

        std::vector<double> zf_vector(N_ + 1, zf_);
        std::vector<double> thetaf_vector(N_ + 1, thetaf_);

        // z_vector and theta_vector from x
        std::vector<double> z_vector(x, x + (N_ + 1)); // Extract the first N+1 elements for z
        std::vector<double> theta_vector(x + 2 * (N_ + 1), x + 3 * (N_ + 1)); // Extract the third N+1 elements for theta

        // resulting difference between z and zf as well as theta and thetaf
        std::vector<double> z_diff(N_ + 1);
        std::vector<double> theta_diff(N_ + 1);
        /*
        // Print zf_vector and thetaf_vector
        std::cout << "zf_vector: ";
        for (const auto& val : zf_vector) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        std::cout << "thetaf_vector: ";
        for (const auto& val : thetaf_vector) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        // Print z_vector and theta_vector
        std::cout << "z_vector: ";
        for (const auto& val : z_vector) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        std::cout << "theta_vector: ";
        for (const auto& val : theta_vector) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        */
        vdSub(N_ + 1, z_vector.data(), zf_vector.data(), z_diff.data()); // z_diff = z_vector - zf_vector
        vdSub(N_ + 1, theta_vector.data(), thetaf_vector.data(), theta_diff.data()); // theta_diff = theta_vector - thetaf_vector
        /*
        // Print z_diff and theta_diff
        std::cout << "z_diff: ";
        for (const auto& val : z_diff) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        std::cout << "theta_diff: ";
        for (const auto& val : theta_diff) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        */
        // calculating the squared differences
        vdSqr(N_ + 1, z_diff.data(), z_diff.data()); // z_diff = (z_vector - zf_vector)^2
        vdSqr(N_ + 1, theta_diff.data(), theta_diff.data()); // theta_diff = (theta_vector - thetaf_vector)^2
        /*
        // Print z_diff^2 and theta_diff^2
        std::cout << "z_diff^2: ";
        for (const auto& val : z_diff) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        std::cout << "theta_diff^2: ";
        for (const auto& val : theta_diff) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        */
        // Sum the squared differences using MKL's cblas_dasum
        double z_diff_sum = cblas_dasum(N_ + 1, z_diff.data(), 1);
        double theta_diff_sum = cblas_dasum(N_ + 1, theta_diff.data(), 1);

        //std::cout << "z_diff_sum: " << z_diff_sum << std::endl;
        //std::cout << "theta_diff_sum: " << theta_diff_sum << std::endl;

        // Compute the weighted sums
        double z_term = w1 * z_diff_sum;
        double theta_term = w2 * theta_diff_sum;

        // Sum the weighted terms to get the objective value
        obj_value = z_term + theta_term;
        //std::cout << "z_term: " << z_term << std::endl;
        //std::cout << "theta_term: " << theta_term << std::endl;
        //std::cout << "objective value = " << obj_value << std::endl;


        return true;
    }



    ///*
    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) {
        // Get the differentiation matrix
        Bebot Bebot(N_, tf_);
        Bebot.calculate();
        const auto& Dm = Bebot.getDifferentiationMatrix();
        // Printing Dm matrix
        //std::cout << "Dm matrix:\n";
        //for (const auto& row : Dm) {
        //    for (const auto& element : row) {
        //        std::cout << element << " ";
        //    }
        //    std::cout << std::endl;
        //}
        // Extracting vectors from x
        std::vector<double> z_vector(x, x + (N_ + 1));
        std::vector<double> w_vector(x + (N_ + 1), x + 2 * (N_ + 1));
        std::vector<double> theta_vector(x + 2 * (N_ + 1), x + 3 * (N_ + 1));
        std::vector<double> q_vector(x + 3 * (N_ + 1), x + 4 * (N_ + 1));
        std::vector<double> delta_v_vector(x + 4 * (N_ + 1), x + 5 * (N_ + 1));
        std::vector<double> delta_s_vector(x + 5 * (N_ + 1), x + 6 * (N_ + 1));
        std::vector<double> delta_m_vector(x + 6 * (N_ + 1), x + 7 * (N_ + 1));
        
        //std::cout << "z_vector: ";
        //for (const auto& val : z_vector) {
        //    std::cout << val << " ";
        //}
        //std::cout << "w_vector: ";
        //for (const auto& val : w_vector) {
        //    std::cout << val << " ";
        //}
        //std::cout << "theta_vector: ";
        //for (const auto& val : theta_vector) {
        //    std::cout << val << " ";
        //}
        //std::cout << "q_vector: ";
        //for (const auto& val : q_vector) {
        //    std::cout << val << " ";
        //}
        //std::cout << "delta_v_vector: ";
        //for (const auto& val : delta_v_vector) {
        //    std::cout << val << " ";
        //}
        //std::cout << "delta_s_vector: ";
        //for (const auto& val : delta_s_vector) {
        //    std::cout << val << " ";
        //}
        //std::cout << "delta_m_vector: ";
        //for (const auto& val : delta_m_vector) {
        //    std::cout << val << " ";
        //}
        
        // Calculate dynamics using differentiation matrix Dm
        std::vector<double> dyn1(N_ + 1);
        std::vector<double> dyn2(N_ + 1);
        std::vector<double> dyn3(N_ + 1);
        std::vector<double> dyn4(N_ + 1);

        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, z_vector.data(), 1, 0.0, dyn1.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, w_vector.data(), 1, 0.0, dyn2.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, theta_vector.data(), 1, 0.0, dyn3.data(), 1);
        cblas_dgemv(CblasColMajor, CblasTrans, N_ + 1, N_ + 1, 1.0, Dm.data(), N_ + 1, q_vector.data(), 1, 0.0, dyn4.data(), 1);
        /*
        std::cout << "dyn1: ";
        for (const auto& val : dyn1) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        std::cout << "dyn2: ";
        for (const auto& val : dyn2) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        std::cout << "dyn3: ";
        for (const auto& val : dyn3) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        
        std::cout << "dyn4: ";
        for (const auto& val : dyn4) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        */
        // Form the X1_matrix and U_matrix
        std::vector<double> X1_matrix_flat(4 * (N_ + 1));
        std::vector<double> U_matrix_flat(3 * (N_ + 1));

        // Fill X1_matrix_flat
        for (Index i = 0; i < (N_ + 1); ++i) {
            X1_matrix_flat[i] = z_vector[i];
            X1_matrix_flat[(N_ + 1) + i] = w_vector[i];
            X1_matrix_flat[2 * (N_ + 1) + i] = theta_vector[i];
            X1_matrix_flat[3 * (N_ + 1) + i] = q_vector[i];
        }

        // Fill U_matrix_flat
        for (Index i = 0; i < (N_ + 1); ++i) {
            U_matrix_flat[i] = delta_v_vector[i];
            U_matrix_flat[(N_ + 1) + i] = delta_s_vector[i];
            U_matrix_flat[2 * (N_ + 1) + i] = delta_m_vector[i];
        }

        // Print X1_matrix_flat
        //std::cout << "X1_matrix_flat: ";
        //for (const auto& val : X1_matrix_flat) {
        //    std::cout << val << " ";
        //}
        //std::cout << std::endl;

        // Print U_matrix_flat
        //std::cout << "U_matrix_flat: ";
        //for (const auto& val : U_matrix_flat) {
        //    std::cout << val << " ";
        //}
        //std::cout << std::endl;

        // Perform matrix multiplication for X2_matrix using MKL
        std::vector<double> X2_matrix_flat(4 * (N_ + 1), 0.0);

        // Form A matrix from A_
        double A[4][4];
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                A[i][j] = A_[i * 4 + j];
            }
        }

        // Form B matrix from B_
        double B[4][3];
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 3; ++j) {
                B[i][j] = B_[i * 3 + j];
            }
        }

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4, (N_ + 1), 4, 1.0, &A[0][0], 4, X1_matrix_flat.data(), (N_ + 1), 0.0, X2_matrix_flat.data(), (N_ + 1));
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4, (N_ + 1), 3, 1.0, &B[0][0], 3, U_matrix_flat.data(), (N_ + 1), 1.0, X2_matrix_flat.data(), (N_ + 1));

        // Print X2_matrix_flat
        //std::cout << "X2_matrix_flat: ";
        //for (const auto& val : X2_matrix_flat) {
        //    std::cout << val << " ";
        //}
        //std::cout << std::endl;

        // Extract z_x2, w_x2, theta_x2, q_x2 from X2_matrix_flat
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
        /*
        std::cout << "z_x2: ";
        for (const auto& val : z_x2) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        std::cout << "w_x2: ";
        for (const auto& val : w_x2) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        std::cout << "theta_x2: ";
        for (const auto& val : theta_x2) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        std::cout << "q_x2: ";
        for (const auto& val : q_x2) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        */
        // Calculate g1, g2, g3, g4
        std::vector<double> g1(N_ + 1);
        std::vector<double> g2(N_ + 1);
        std::vector<double> g3(N_ + 1);
        std::vector<double> g4(N_ + 1);

        vdSub(N_ + 1, dyn1.data(), z_x2.data(), g1.data());
        vdSub(N_ + 1, dyn2.data(), w_x2.data(), g2.data());
        vdSub(N_ + 1, dyn3.data(), theta_x2.data(), g3.data());
        vdSub(N_ + 1, dyn4.data(), q_x2.data(), g4.data());

        // Fill the g vector
        for (Index i = 0; i < (N_ + 1); ++i) {
            g[i] = g1[i];
            g[(N_ + 1) + i] = g2[i];
            g[2 * (N_ + 1) + i] = g3[i];
            g[3 * (N_ + 1) + i] = g4[i];
        }

        // Print g vector
        //std::cout << "g vector: ";
        //for (Index i = 0; i < 4 * (N_ + 1); ++i) {
        //    std::cout << g[i] << " ";
        //}
        //std::cout << std::endl;

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

    ) { // Resize solution vectors
        solution_x_.resize(4 * (N_ + 1));
        solution_u_.resize(3 * (N_ + 1));
        
        // Extract solution_x_ and solution_u_ from x
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
        
                // Create BernsteinPoly results for each variable
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
   
        // Flatten results for saving
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
    // Getter for the solution
    const std::vector<Number>& get_solution_x() const { return solution_x_; }

    // Getter for the final objective value
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
    std::vector<double> A_;
    std::vector<double> B_;


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

int main() {
    int N = 30;
    double tf = 10.0;
    double delta_v_max = 30.0;
    double delta_v_min = -30.0;
    double delta_s_max = 30.0;
    double delta_s_min = -30.0;
    double delta_m_max = 4000.0;
    double delta_m_min = -4000.0;
    double zmax = 0.0;
    double zmin = -100.0;
    double wmax = 5.0;
    double wmin = -5.0;
    double thetamax = 0.5;
    double thetamin = -0.5;
    double qmax = 5.0;
    double qmin = -5.0;
    double z0 = -20.0;
    double w0 = 0.1;
    double theta0 = 0.01;
    double q0 = 0.01;
    double delta_v0 = -20;
    double delta_s0 = -5.0;
    double delta_m0 = -300.00;

    double zf = -25.0;
    double thetaf = 0.0;
    double speed = 2.72;

    std::vector<double> A_flat;
    std::vector<double> B_flat;

    getStateSpaceMatrices(speed, A_flat, B_flat);

    SmartPtr<TNLP> pointSetProblem = new PointSetProblem(N, tf, delta_v_max, delta_v_min, delta_s_max, delta_s_min, delta_m_max, delta_m_min, zmax, zmin, wmax, wmin, thetamax, thetamin, qmax, qmin, z0, w0, theta0, q0, delta_v0, delta_s0, delta_m0, zf, thetaf, A_flat, B_flat);
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
    app->Options()->SetNumericValue("tol", 1e-3); // Change to my desired tolerance



    app->RethrowNonIpoptException(true);
    ApplicationReturnStatus status = app->Initialize();
    if (status != Solve_Succeeded) {
        std::cout << "IPOPT initialization failed!" << std::endl;
        return -1;
    }

    status = app->OptimizeTNLP(pointSetProblem);

    if (status == Solve_Succeeded || status == Solved_To_Acceptable_Level) {
        // Retrieve the optimal solution and objective value from the problem
        const auto& solution_x = static_cast<PointSetProblem*>(GetRawPtr(pointSetProblem))->get_solution_x();
        Number final_obj_value = static_cast<PointSetProblem*>(GetRawPtr(pointSetProblem))->get_final_obj_value();
    
        std::cout << "Optimal Solution (x): ";
        for (Index i = 0; i < solution_x.size(); i++) {
            std::cout << solution_x[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "Optimal Objective Value: " << final_obj_value << std::endl;
        
        // Retrieve the Bernstein polynomial result from the problem
        const auto& bernsteinpoly_result = static_cast<PointSetProblem*>(GetRawPtr(pointSetProblem))->get_bernsteinpoly_result();
        //std::cout << "Bernstein Polynomial Result (xN): " << std::endl;
        //for (const auto& row : bernsteinpoly_result) {
        //    for (double value : row) {
        //        std::cout << value << " ";
        //    }
        //    std::cout << std::endl; // New line for each row
        //}
    } else {
        std::cout << "IPOPT optimization failed with status " << status << std::endl;
    }

    return 0;
}

// Compile and run instructions
// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp_v2/examples/bebot/ma_57/example_mpc_auv_v1$ g++ -o bebot_example_mpc_auv_v1 ~/dev/optimization/BeBOT_cpp_v2/examples/bebot/ma_57/example_mpc_auv_v1/bebot_example_mpc_auv_v1.cpp ~/dev/optimization/BeBOT_cpp_v2/examples/bebot/ma_57/example_mpc_auv_v1/state_space_matrices.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bebot.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinpoly.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteindifferentialmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/bernsteinmatrix_a2b.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/degelevmatrix.cpp ~/dev/optimization/BeBOT_cpp_v2/bebot/nchoosek_mod.cpp -I~/dev/optimization/BeBOT_cpp_v2/include -I./Ipopt/src/ -L./Ipopt/src/.libs -lipopt -L/opt/intel/oneapi/mkl/latest/lib/intel64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -ldl -lm -lpthread -lstdc++
// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp_v2/examples/bebot/ma_57/example_mpc_auv_v1$ export LD_LIBRARY_PATH=/usr/local/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
// vpetrov@lnx-me002:~/dev/optimization/BeBOT_cpp_v2/examples/bebot/ma_57/example_mpc_auv_v1$ ./bebot_example_mpc_auv_v1
