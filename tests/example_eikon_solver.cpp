#include "eikon_solver_2d.hpp"
#include "hdf5_utils.hpp"
#include <iostream>
#include <chrono>

int main() {
    // Example usage of the C++ Eikonal solver
    
    // Define grid - use smaller grid for testing
    int nx = 101, ny = 101;  // 减小网格尺寸
    double x_min = -2.0, x_max = 2.0;  // 减小计算域
    double y_min = -2.0, y_max = 2.0;
    
    // Create coordinate vectors
    Eigen::VectorXd xx = Eigen::VectorXd::LinSpaced(nx, x_min, x_max);
    Eigen::VectorXd yy = Eigen::VectorXd::LinSpaced(ny, y_min, y_max);
    
    // Initialize coefficient matrices
    Eigen::MatrixXd a(nx, ny);
    Eigen::MatrixXd b(nx, ny);
    Eigen::MatrixXd c(nx, ny);
    Eigen::MatrixXd fun(nx, ny);  // slowness field
    Eigen::MatrixXd u(nx, ny);    // reference solution (not used in computation)
    
    // Example: homogeneous medium with anisotropy  
    double a_val = 1.0, b_val = 1.0, c_val = 0.0;  // 确保是各向同性
    double slowness = 1.0;  // 使用更简单的slowness = 1

    a.setConstant(a_val);
    b.setConstant(b_val);
    c.setConstant(c_val);
    fun.setConstant(slowness);
    u.setZero();  // Initialize reference solution
    
    // Source location
    double x0 = 0.0, y0 = 0.0;
    
    // Output matrix
    Eigen::MatrixXd T(nx, ny);
    
    // std::cout << "Starting Eikonal equation solver..." << std::endl;
    // std::cout << "Grid size: " << nx << " x " << ny << std::endl;
    // std::cout << "Domain: [" << x_min << ", " << x_max << "] x [" << y_min << ", " << y_max << "]" << std::endl;
    // std::cout << "Source location: (" << x0 << ", " << y0 << ")" << std::endl;
    // std::cout << "Coefficients: a=" << a_val << ", b=" << b_val << ", c=" << c_val << std::endl;
    // std::cout << "Slowness: " << slowness << " (velocity = " << 1.0/slowness << ")" << std::endl;
    
    // Time the solver
    auto start = std::chrono::high_resolution_clock::now();
    
    // Solve the Eikonal equation
    EikonSolver2D::FSM_WENO3_PS_2d(xx, yy, a, b, c, fun, x0, y0, T);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "Solver completed in " << duration.count() << " ms" << std::endl;
    
    // 添加调试信息
    std::cout << "\nDebugging solver output:" << std::endl;
    std::cout << "Expected max travel time for this domain: " << slowness * std::sqrt(4.0 + 4.0) << std::endl;
    std::cout << "T matrix near source:" << std::endl;
    int mid_x = nx / 2, mid_y = ny / 2;
    for (int i = mid_x - 2; i <= mid_x + 2; ++i) {
        for (int j = mid_y - 2; j <= mid_y + 2; ++j) {
            if (i >= 0 && i < nx && j >= 0 && j < ny) {
                std::cout << "T(" << i << "," << j << ") = " << T(i, j) << "  ";
            }
        }
        std::cout << std::endl;
    }
    
    // 检查几个边界点的值
    std::cout << "\nBoundary points analysis:" << std::endl;
    std::cout << "T(0,0) = " << T(0, 0) << " (corner)" << std::endl;
    std::cout << "T(" << nx-1 << "," << ny-1 << ") = " << T(nx-1, ny-1) << " (opposite corner)" << std::endl;
    std::cout << "T(0," << mid_y << ") = " << T(0, mid_y) << " (left edge center)" << std::endl;
    std::cout << "T(" << nx-1 << "," << mid_y << ") = " << T(nx-1, mid_y) << " (right edge center)" << std::endl;
    
    // Compute analytical solution for comparison
    std::cout << "\nComputing analytical solution for comparison..." << std::endl;
    Eigen::MatrixXd T_analytical(nx, ny);
    
    // Calculate grid spacing
    double dx = (x_max - x_min) / (nx - 1);
    double dy = (y_max - y_min) / (ny - 1);
    
    std::cout << "Grid spacing: dx = " << dx << ", dy = " << dy << std::endl;
    std::cout << "Slowness = " << slowness << " (velocity = " << 1.0/slowness << ")" << std::endl;
    
    // For isotropic medium (c=0), the analytical solution is the Euclidean distance
    // multiplied by slowness: T = slowness * sqrt((x-x0)^2 + (y-y0)^2)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double dx_coord = xx(i) - x0;
            double dy_coord = yy(j) - y0;
            T_analytical(i, j) = slowness * std::sqrt(dx_coord*dx_coord + dy_coord*dy_coord);
        }
    }
  
    // Print travel times at a few sample points
    std::cout << "\nSample travel times comparison:" << std::endl;
    
    // Check points at various distances from source
    std::vector<std::pair<int, int>> test_points = {
        {mid_x, mid_y},           // source point
        {mid_x + 1, mid_y},       // 1 grid point away
        {mid_x + 2, mid_y},       // 2 grid points away
        {mid_x + 5, mid_y},       // 5 grid points away
        {mid_x + 10, mid_y},      // 10 grid points away
        {mid_x, mid_y + 1},       // 1 grid point away in y
        {mid_x + 1, mid_y + 1}    // diagonal
    };
    
    for (auto& point : test_points) {
        int i = point.first, j = point.second;
        if (i >= 0 && i < nx && j >= 0 && j < ny) {
            double distance = std::sqrt(std::pow(xx(i) - x0, 2) + std::pow(yy(j) - y0, 2));
            std::cout << "At (" << i << ", " << j << ") - distance " << distance << ":" << std::endl;
            std::cout << "  Numerical: " << T(i, j) << std::endl;
            std::cout << "  Analytical: " << T_analytical(i, j) << std::endl;
            std::cout << "  Error: " << std::abs(T(i, j) - T_analytical(i, j)) << std::endl;
            std::cout << "  Rel Error: " << (std::abs(T(i, j) - T_analytical(i, j)) / std::max(T_analytical(i, j), 1e-12)) * 100 << "%" << std::endl;
            std::cout << std::endl;
        }
    }
    
    // Save results to HDF5 using the dedicated utility class
    std::cout << "\nSaving results to HDF5..." << std::endl;
    
    // Save only travel time matrix with coordinates
    HDF5Utils::saveTravelTimeResults(xx, yy, T, "travel_times.h5");
    
    // Save complete results including coefficients and source location
    HDF5Utils::saveCompleteResults(xx, yy, T, "complete_results.h5", &a, &b, &c, &fun, &x0, &y0);
    
    // Alternative: save just the T matrix if that's all you need
    HDF5Utils::saveMatrix(T, "T_matrix_only.h5", "travel_time");
    
    std::cout << "\nExample completed successfully!" << std::endl;
    
    return 0;
}
