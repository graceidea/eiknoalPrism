#include "eikon_solver_2d.hpp"

void EikonSolver2D::FSM_WENO3_PS_2d(
    const Eigen::VectorXd& xx,
    const Eigen::VectorXd& yy,
    const Eigen::MatrixXd& a,
    const Eigen::MatrixXd& b,
    const Eigen::MatrixXd& c,
    const Eigen::MatrixXd& fun,
    double x0,
    double y0,
    Eigen::MatrixXd& T,
    bool use_cuda
) {
    int nx = static_cast<int>(xx.size());
    int ny = static_cast<int>(yy.size());
    
    // Initialize output matrix
    T.resize(nx, ny);
    
    // Create mesh spacing
    double dx = xx(1) - xx(0);
    double dy = yy(1) - yy(0);
    
    // Parameter discretization at the seismic source
    int idx0 = static_cast<int>(std::floor((x0 - xx(0)) / dx));
    int idy0 = static_cast<int>(std::floor((y0 - yy(0)) / dy));
    
    // Ensure indices are within bounds
    idx0 = std::max(0, std::min(idx0, nx - 2));
    idy0 = std::max(0, std::min(idy0, ny - 2));
    
    double r1 = std::min(1.0, (x0 - xx(idx0)) / dx);
    double r2 = std::min(1.0, (y0 - yy(idy0)) / dy);
    
    // Bilinear interpolation at source location
    double a0 = (1 - r1) * (1 - r2) * a(idx0, idy0) + 
                (1 - r1) * r2 * a(idx0, idy0 + 1) +
                r1 * (1 - r2) * a(idx0 + 1, idy0) + 
                r1 * r2 * a(idx0 + 1, idy0 + 1);
    
    double b0 = (1 - r1) * (1 - r2) * b(idx0, idy0) + 
                (1 - r1) * r2 * b(idx0, idy0 + 1) +
                r1 * (1 - r2) * b(idx0 + 1, idy0) + 
                r1 * r2 * b(idx0 + 1, idy0 + 1);
    
    double c0 = (1 - r1) * (1 - r2) * c(idx0, idy0) + 
                (1 - r1) * r2 * c(idx0, idy0 + 1) +
                r1 * (1 - r2) * c(idx0 + 1, idy0) + 
                r1 * r2 * c(idx0 + 1, idy0 + 1);
    
    double fun0 = (1 - r1) * (1 - r2) * fun(idx0, idy0) + 
                  (1 - r1) * r2 * fun(idx0, idy0 + 1) +
                  r1 * (1 - r2) * fun(idx0 + 1, idy0) + 
                  r1 * r2 * fun(idx0 + 1, idy0 + 1);
    
    // Initialize matrices
    Eigen::MatrixXd T0(nx, ny), T0x(nx, ny), T0y(nx, ny);
    Eigen::MatrixXd tau(nx, ny);
    Eigen::MatrixXi ischange(nx, ny);
    
    // Compute T0 and its derivatives
    computeT0AndDerivatives(xx, yy, x0, y0, a0, b0, c0, fun0, T0, T0x, T0y);
    
    // Initialize tau and ischange arrays
    initializeTauAndChange(xx, yy, x0, y0, dx, dy, tau, ischange);
    
    // Main iteration loop
    Eigen::MatrixXd tau_old(nx, ny);
    double L1_dif, Linf_dif, L1_err, Linf_err;
    
    // Pre-compute vectorized matrices for efficiency
    Eigen::MatrixXd sigx_matrix = a.array().sqrt();
    Eigen::MatrixXd sigy_matrix = b.array().sqrt();
    Eigen::MatrixXd inv_dx_dy = (sigx_matrix / dx + sigy_matrix / dy).cwiseInverse();


#ifdef CUDA_ENABLED
    // Initialize GPU if requested and available
    if (use_cuda) {
        try {
            initializeGPU(nx, ny);
            copyDataToGPU(a, b, c, fun, T0, T0x, T0y, ischange, tau, dx, dy);
        } catch (const std::exception& e) {
            std::cerr << "GPU initialization failed: " << e.what() << std::endl;
            std::cerr << "Falling back to CPU computation" << std::endl;
            use_cuda = false;
        }
        // GPU-accelerated sweep
        sweepAllDirectionsGPU(T, dx, dy, nx, ny);
    }
#else
    use_cuda = false;
#endif
    if (!use_cuda) {
        for (int iter = 0; iter < MaxIter; ++iter) {
            tau_old = tau;
            // Four-direction sweeping
            for (int xdirec = -1; xdirec <= 1; xdirec += 2) {
                for (int ydirec = -1; ydirec <= 1; ydirec += 2) {

                    // Determine sweep order based on direction
                    int x_start, x_end, y_start, y_end;
                    selectDirection(xdirec, ydirec, nx, ny,
                                    x_start, x_end, y_start, y_end);

                    // Choose GPU or CPU computation
                    // CPU sweep (original implementation)
                    // Sweep in the specified direction
                    for (int iix = x_start; iix != x_end; iix -= xdirec) {
                        for (int iiy = y_start; iiy != y_end; iiy -= ydirec) {

                            if (ischange(iix, iiy) == 1) {
                                // Use pre-computed values (vectorized lookup)
                                double sigx = sigx_matrix(iix, iiy);
                                double sigy = sigy_matrix(iix, iiy);
                                double coe = inv_dx_dy(iix, iiy);
                                
                                // Compute WENO derivatives
                                double px1, px2, py1, py2;
                                computeWENODerivatives(tau, iix, iiy, nx, ny, dx, dy, px1, px2, py1, py2);
                                
                                // Optimized Hamiltonian computation using pre-computed terms
                                double px_avg = (px1 + px2) * 0.5;
                                double py_avg = (py1 + py2) * 0.5;
                                
                                double Htau = std::sqrt(
                                    a(iix, iiy) * px_avg * px_avg + 
                                    b(iix, iiy) * py_avg * py_avg -
                                    2.0 * c(iix, iiy) * px_avg * py_avg +
                                    2.0 * (a(iix, iiy) * T0x(iix, iiy) - c(iix, iiy) * T0y(iix, iiy)) * px_avg +
                                    2.0 * (b(iix, iiy) * T0y(iix, iiy) - c(iix, iiy) * T0x(iix, iiy)) * py_avg +
                                    a(iix, iiy) * T0x(iix, iiy) * T0x(iix, iiy) + 
                                    b(iix, iiy) * T0y(iix, iiy) * T0y(iix, iiy) -
                                    2.0 * c(iix, iiy) * T0x(iix, iiy) * T0y(iix, iiy)
                                );
                                
                                // Update tau using correct Fast Sweeping formula
                                // Remove the incorrect addition of current tau value
                                double tpT = coe * (fun(iix, iiy) - Htau) + 
                                            coe * (sigx * (px2 - px1) / 2.0 + sigy * (py2 - py1) / 2.0) + 
                                            tau(iix, iiy);
                                
                                // Only update if the new value is smaller (Fast Sweeping Method principle)
                                tau(iix, iiy) = tpT;
                            }
                        }
                    }
                    // Update boundary conditions
                    updateBoundaryConditions(tau, nx, ny);
                }
                
            }
            // Compute convergence metrics
            computeConvergenceMetrics(tau, tau_old, T0, dx, dy, nx, ny, 
                                 L1_dif, Linf_dif, L1_err, Linf_err);
            // Check convergence
            if (std::abs(L1_dif) < tol && Linf_dif < tol) {
                break;
            }
        }
        // Final result: T = tau + T0
        T = tau + T0;
    }

}

void EikonSolver2D::computeT0AndDerivatives(
    const Eigen::VectorXd& xx,
    const Eigen::VectorXd& yy,
    double x0, double y0,
    double a0, double b0, double c0, double fun0,
    Eigen::MatrixXd& T0,
    Eigen::MatrixXd& T0x,
    Eigen::MatrixXd& T0y
) {
    int nx = static_cast<int>(xx.size());
    int ny = static_cast<int>(yy.size());
    
    double det = a0 * b0 - c0 * c0;
    double inv_det = 1.0 / det;
    
    // Create grid matrices using broadcasting == meshgrid
    Eigen::MatrixXd X = xx.replicate(1, ny);
    Eigen::MatrixXd Y = yy.transpose().replicate(nx, 1);
    
    // Compute displacement from source (vectorized)
    Eigen::MatrixXd dx_grid = X.array() - x0;
    Eigen::MatrixXd dy_grid = Y.array() - y0;
    
    // Compute quadratic form for all points (vectorized)
    Eigen::MatrixXd quadratic_form = 
        b0 * inv_det * dx_grid.array().square() + 
        a0 * inv_det * dy_grid.array().square() + 
        2.0 * c0 * inv_det * dx_grid.cwiseProduct(dy_grid).array();
    
    // Compute T0 (vectorized)
    T0 = fun0 * quadratic_form.array().sqrt();
    
    // Create mask for zero T0 values
    Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> zero_mask = 
        (T0.array() == 0.0);
    
    // Compute derivatives for non-zero T0 values (vectorized)
    Eigen::MatrixXd T0x_nonzero = 
        fun0 * fun0 * (b0 * inv_det * dx_grid.array() + c0 * inv_det * dy_grid.array()) / T0.array();
    Eigen::MatrixXd T0y_nonzero = 
        fun0 * fun0 * (a0 * inv_det * dy_grid.array() + c0 * inv_det * dx_grid.array()) / T0.array();
    
    // Apply conditional assignment using select (vectorized)
    T0x = zero_mask.select(0.0, T0x_nonzero);
    T0y = zero_mask.select(0.0, T0y_nonzero);
}

void EikonSolver2D::initializeTauAndChange(
    const Eigen::VectorXd& xx,
    const Eigen::VectorXd& yy,
    double x0, double y0,
    double dx, double dy,
    Eigen::MatrixXd& tau,
    Eigen::MatrixXi& ischange
) {
    int nx = static_cast<int>(xx.size());
    int ny = static_cast<int>(yy.size());
    
    // Create grid matrices using broadcasting
    Eigen::MatrixXd X = xx.replicate(1, ny);
    Eigen::MatrixXd Y = yy.transpose().replicate(nx, 1);
    
    // Compute normalized grid distances
    Eigen::MatrixXd dx_grid = (X.array() - x0) / dx;
    Eigen::MatrixXd dy_grid = (Y.array() - y0) / dy;
    
    // Create mask for source region (where abs(dx_grid) <= 0.9 && abs(dy_grid) <= 0.9)
    Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> source_mask = 
        (dx_grid.array().abs() <= 0.9) && (dy_grid.array().abs() <= 0.9);
    
    // Initialize tau and ischange using vectorized operations
    tau = Eigen::MatrixXd::Constant(nx, ny, TAU_INF);
    ischange = Eigen::MatrixXi::Ones(nx, ny);
    
    // Set source region values using vectorized operations
    tau = source_mask.select(0.0, tau);
    ischange = source_mask.select(0, ischange);
    
    // Check if source is on boundary (vectorized boundary detection)
    Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> boundary_mask = 
        (Eigen::ArrayXi::LinSpaced(nx, 0, nx-1).replicate(1, ny) == 0) ||
        (Eigen::ArrayXi::LinSpaced(nx, 0, nx-1).replicate(1, ny) == nx-1) ||
        (Eigen::ArrayXi::LinSpaced(ny, 0, ny-1).transpose().replicate(nx, 1) == 0) ||
        (Eigen::ArrayXi::LinSpaced(ny, 0, ny-1).transpose().replicate(nx, 1) == ny-1);
    
    if ((source_mask && boundary_mask).any()) {
        std::cerr << "Warning: source on the boundary, mesh error" << std::endl;
    }
}

void EikonSolver2D::computeWENODerivatives(
    const Eigen::MatrixXd& tau,
    int iix, int iiy, int nx, int ny,
    double dx, double dy,
    double& px1, double& px2,
    double& py1, double& py2
) {
    // X-direction derivatives
    if (iix == 1) {
        px1 = (tau(iix, iiy) - tau(iix - 1, iiy)) / dx;
        
        double numerator = eps + std::pow(tau(iix, iiy) - 2 * tau(iix + 1, iiy) + tau(iix + 2, iiy), 2);
        double denominator = eps + std::pow(tau(iix - 1, iiy) - 2 * tau(iix, iiy) + tau(iix + 1, iiy), 2);
        double wx2 = 1.0 / (1.0 + 2.0 * std::pow(numerator / denominator, 2));
        
        px2 = (1 - wx2) * (tau(iix + 1, iiy) - tau(iix - 1, iiy)) / (2.0 * dx) +
              wx2 * (-3 * tau(iix, iiy) + 4 * tau(iix + 1, iiy) - tau(iix + 2, iiy)) / (2.0 * dx);
    }
    else if (iix == nx - 2) {
        double numerator = eps + std::pow(tau(iix, iiy) - 2 * tau(iix - 1, iiy) + tau(iix - 2, iiy), 2);
        double denominator = eps + std::pow(tau(iix + 1, iiy) - 2 * tau(iix, iiy) + tau(iix - 1, iiy), 2);
        double wx1 = 1.0 / (1.0 + 2.0 * std::pow(numerator / denominator, 2));
        
        px1 = (1 - wx1) * (tau(iix + 1, iiy) - tau(iix - 1, iiy)) / (2.0 * dx) +
              wx1 * (3 * tau(iix, iiy) - 4 * tau(iix - 1, iiy) + tau(iix - 2, iiy)) / (2.0 * dx);
        px2 = (tau(iix + 1, iiy) - tau(iix, iiy)) / dx;
    }
    else {
        // Interior points
        double numerator1 = eps + std::pow(tau(iix, iiy) - 2 * tau(iix - 1, iiy) + tau(iix - 2, iiy), 2);
        double denominator1 = eps + std::pow(tau(iix + 1, iiy) - 2 * tau(iix, iiy) + tau(iix - 1, iiy), 2);
        double wx1 = 1.0 / (1.0 + 2.0 * std::pow(numerator1 / denominator1, 2));
        
        px1 = (1.0 - wx1) * (tau(iix + 1, iiy) - tau(iix - 1, iiy)) / (2.0 * dx) +
              wx1 * (3 * tau(iix, iiy) - 4 * tau(iix - 1, iiy) + tau(iix - 2, iiy)) / (2.0 * dx);
        
        double numerator2 = eps + std::pow(tau(iix, iiy) - 2 * tau(iix + 1, iiy) + tau(iix + 2, iiy), 2);
        double denominator2 = eps + std::pow(tau(iix - 1, iiy) - 2 * tau(iix, iiy) + tau(iix + 1, iiy), 2);
        double wx2 = 1.0 / (1.0 + 2.0 * std::pow(numerator2 / denominator2, 2));
        
        px2 = (1.0 - wx2) * (tau(iix + 1, iiy) - tau(iix - 1, iiy)) / (2.0 * dx) +
              wx2 * (-3 * tau(iix, iiy) + 4 * tau(iix + 1, iiy) - tau(iix + 2, iiy)) / (2.0 * dx);
    }
    
    // Y-direction derivatives
    if (iiy == 1) {
        py1 = (tau(iix, iiy) - tau(iix, iiy - 1)) / dy;
        
        double numerator = eps + std::pow(tau(iix, iiy) - 2 * tau(iix, iiy + 1) + tau(iix, iiy + 2), 2);
        double denominator = eps + std::pow(tau(iix, iiy - 1) - 2 * tau(iix, iiy) + tau(iix, iiy + 1), 2);
        double wy2 = 1.0 / (1.0 + 2.0 * std::pow(numerator / denominator, 2));
        
        py2 = (1 - wy2) * (tau(iix, iiy + 1) - tau(iix, iiy - 1)) / (2.0 * dy) +
              wy2 * (-3 * tau(iix, iiy) + 4 * tau(iix, iiy + 1) - tau(iix, iiy + 2)) / (2.0 * dy);
    }
    else if (iiy == ny - 2) {
        double numerator = eps + std::pow(tau(iix, iiy) - 2 * tau(iix, iiy - 1) + tau(iix, iiy - 2), 2);
        double denominator = eps + std::pow(tau(iix, iiy + 1) - 2 * tau(iix, iiy) + tau(iix, iiy - 1), 2);
        double wy1 = 1.0 / (1.0 + 2.0 * std::pow(numerator / denominator, 2));
        
        py1 = (1 - wy1) * (tau(iix, iiy + 1) - tau(iix, iiy - 1)) / (2.0 * dy) +
              wy1 * (3 * tau(iix, iiy) - 4 * tau(iix, iiy - 1) + tau(iix, iiy - 2)) / (2.0 * dy);
        py2 = (tau(iix, iiy + 1) - tau(iix, iiy)) / dy;
    }
    else {
        // Interior points
        double numerator1 = eps + std::pow(tau(iix, iiy) - 2 * tau(iix, iiy - 1) + tau(iix, iiy - 2), 2);
        double denominator1 = eps + std::pow(tau(iix, iiy + 1) - 2 * tau(iix, iiy) + tau(iix, iiy - 1), 2);
        double wy1 = 1.0 / (1.0 + 2.0 * std::pow(numerator1 / denominator1, 2));
        
        py1 = (1 - wy1) * (tau(iix, iiy + 1) - tau(iix, iiy - 1)) / (2.0 * dy) +
              wy1 * (3 * tau(iix, iiy) - 4 * tau(iix, iiy - 1) + tau(iix, iiy - 2)) / (2.0 * dy);
        
        double numerator2 = eps + std::pow(tau(iix, iiy) - 2 * tau(iix, iiy + 1) + tau(iix, iiy + 2), 2);
        double denominator2 = eps + std::pow(tau(iix, iiy - 1) - 2 * tau(iix, iiy) + tau(iix, iiy + 1), 2);
        double wy2 = 1.0 / (1.0 + 2.0 * std::pow(numerator2 / denominator2, 2));
        
        py2 = (1 - wy2) * (tau(iix, iiy + 1) - tau(iix, iiy - 1)) / (2.0 * dy) +
              wy2 * (-3 * tau(iix, iiy) + 4 * tau(iix, iiy + 1) - tau(iix, iiy + 2)) / (2.0 * dy);
    }
}

void EikonSolver2D::updateBoundaryConditions(
    Eigen::MatrixXd& tau,
    int nx, int ny
) {
    // Update boundary conditions using extrapolation (vectorized)
    // Left and right boundaries
    tau.col(0) = (2.0 * tau.col(1) - tau.col(2)).cwiseMax(tau.col(2));
    tau.col(nx-1) = (2.0 * tau.col(nx-2) - tau.col(nx-3)).cwiseMax(tau.col(nx-3));
    
    // Top and bottom boundaries  
    tau.row(0) = (2.0 * tau.row(1) - tau.row(2)).cwiseMax(tau.row(2));
    tau.row(ny-1) = (2.0 * tau.row(ny-2) - tau.row(ny-3)).cwiseMax(tau.row(ny-3));
}

void EikonSolver2D::computeConvergenceMetrics(
    const Eigen::MatrixXd& tau,
    const Eigen::MatrixXd& tau_old,
    const Eigen::MatrixXd& T0,
    double dx, double dy,
    int nx, int ny,
    double& L1_dif, double& Linf_dif,
    double& L1_err, double& Linf_err
) {
    L1_dif = 0.0;
    Linf_dif = 0.0;
    L1_err = 0.0;
    Linf_err = 0.0;
    
    // Compute L1 and L_infinity differences
    L1_dif = (tau - tau_old).cwiseAbs().sum() * dx * dy;
    Linf_dif = (tau - tau_old).cwiseAbs().maxCoeff();
    
    // Compute L1 and L_infinity errors (excluding boundary regions)
    // int count = 0;
    // auto ind = Eigen::seq(4, nx-5);
    // L1_err = (tau(ind, ind) + T0(ind, ind) - ref_t).cwiseAbs().sum() * dx * dy;
    // Linf_err = (tau(ind, ind) + T0(ind, ind) - ref_t).cwiseAbs().maxCoeff();
    
    // if (count > 0) {
    //     L1_err /= count;
    // }
}

void EikonSolver2D::selectDirection(
    const int xdirec,
    const int ydirec,
    int nx, int ny,
    int& x_start, int& x_end, int& y_start, int& y_end
) {
    // Select the sweeping direction
    if (xdirec < 0) {  // Right
        x_start = 1;
        x_end = nx - 1;
    } else {  // Left
        x_start = nx - 2;
        x_end = 2;
    }

    if (ydirec < 0) {  // Down
        y_start = 1;
        y_end = ny - 1;
    } else {  // Up
        y_start = ny - 2;
        y_end = 2;
    }
}