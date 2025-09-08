#ifndef EIKON_SOLVER_2D_HPP
#define EIKON_SOLVER_2D_HPP

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <algorithm>

class EikonSolver2D {
public:
    /**
     * FSM_WENO3_PS_2d - Fast Sweeping Method with WENO3 scheme for 2D Eikonal equation
     * 
     * @param xx: x coordinates vector (nx elements)
     * @param yy: y coordinates vector (ny elements) 
     * @param a: coefficient matrix a (nx x ny)
     * @param b: coefficient matrix b (nx x ny)
     * @param c: coefficient matrix c (nx x ny)
     * @param fun: slowness field (nx x ny)
     * @param x0: source x coordinate
     * @param y0: source y coordinate
     * @param u: reference solution (nx x ny) - not used in computation
     * @param T: output travel time field (nx x ny)
     * @param use_cuda: whether to use CUDA GPU acceleration (default: false)
     */
    static void FSM_WENO3_PS_2d(
        const Eigen::VectorXd& xx,
        const Eigen::VectorXd& yy,
        const Eigen::MatrixXd& a,
        const Eigen::MatrixXd& b,
        const Eigen::MatrixXd& c,
        const Eigen::MatrixXd& fun,
        double x0,
        double y0,
        Eigen::MatrixXd& T,
        bool use_cuda = false
    );

private:
    static constexpr int MaxIter = 1000;
    static constexpr double tol = 1e-4;
    static constexpr double eps = 1e-14;
    static constexpr double TAU_INF = 20.0;
    static constexpr double ref_t = 1.0;

    // Helper functions
    static void selectDirection(
        const int xdirec,
        const int ydirec,
        int nx, int ny,
        int& x_start, int& x_end, int& y_start, int& y_end
    );

    static void computeT0AndDerivatives(
        const Eigen::VectorXd& xx,
        const Eigen::VectorXd& yy,
        double x0, double y0,
        double a0, double b0, double c0, double fun0,
        Eigen::MatrixXd& T0,
        Eigen::MatrixXd& T0x,
        Eigen::MatrixXd& T0y
    );
    
    static void initializeTauAndChange(
        const Eigen::VectorXd& xx,
        const Eigen::VectorXd& yy,
        double x0, double y0,
        double dx, double dy,
        Eigen::MatrixXd& tau,
        Eigen::MatrixXi& ischange
    );
    
    static void computeWENODerivatives(
        const Eigen::MatrixXd& tau,
        int iix, int iiy, int nx, int ny,
        double dx, double dy,
        double& px1, double& px2,
        double& py1, double& py2
    );
    
    static void updateBoundaryConditions(
        Eigen::MatrixXd& tau,
        int nx, int ny
    );
    
    static void computeConvergenceMetrics(
        const Eigen::MatrixXd& tau,
        const Eigen::MatrixXd& tau_old,
        const Eigen::MatrixXd& T0,
        double dx, double dy,
        int nx, int ny,
        double& L1_dif, double& Linf_dif,
        double& L1_err, double& Linf_err
    );

#ifdef CUDA_ENABLED
    // CUDA GPU acceleration functions
    static void initializeGPU(int nx, int ny);
    static void copyDataToGPU(
        const Eigen::MatrixXd& a, const Eigen::MatrixXd& b, 
        const Eigen::MatrixXd& c, const Eigen::MatrixXd& fun,
        const Eigen::MatrixXd& T0, const Eigen::MatrixXd& T0x,
        const Eigen::MatrixXd& T0y, const Eigen::MatrixXi& ischange,
        const Eigen::MatrixXd& tau,
        double dx, double dy
    );
   
    static void computeSimpleConvergenceMetricsGPU(
        int nx, int ny, double dx, double dy, double& L1_dif, double& Linf_dif
    );
   
    static void sweepAllDirectionsGPU(
        Eigen::MatrixXd& T, double dx, double dy, int nx, int ny
    );
    static void cleanupGPU();
    
    // CUDA device memory pointers (static for class-level management)
    static double* d_tau;
    static double* d_tau_old;
    static double* d_T0;
    static double* d_T0x;
    static double* d_T0y;
    static double* d_a;
    static double* d_b;
    static double* d_c;
    static double* d_fun;
    static int* d_ischange;
    static double* d_sigx_matrix;
    static double* d_sigy_matrix;
    static double* d_inv_dx_dy;
    static double* d_H_const;
    static double* d_H_coeff_x;
    static double* d_H_coeff_y;
    static double* d_reduction_buffer;
    static double* d_diff_buffer; 
    static double* d_T;
    static bool gpu_initialized;
#endif
};

#endif // EIKON_SOLVER_2D_HPP
