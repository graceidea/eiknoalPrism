#include "eikon_solver_2d.hpp"
#include "hdf5_utils.hpp"
#include <iostream>
#include <chrono>

// Analytical solution for isotropic case (for validation)
double analytical_isotropic(double x, double y, double x0, double y0, double vel) {
    double dx = x - x0;
    double dy = y - y0;
    return std::sqrt(dx*dx + dy*dy) / vel;
}

int main() {
    try {
        std::cout << "=== Eikonal Solver GPU vs CPU Comparison ===" << std::endl;
        
        // Create test problem (isotropic case for validation)
        int nx = 201;
        int ny = 201;
        double xmin = 0.0, xmax = 10.0;
        double ymin = 0.0, ymax = 10.0;
        
        std::cout << "Grid size: " << nx << " x " << ny << std::endl;
        
        // Create coordinate vectors
        Eigen::VectorXd xx = Eigen::VectorXd::LinSpaced(nx, xmin, xmax);
        Eigen::VectorXd yy = Eigen::VectorXd::LinSpaced(ny, ymin, ymax);
        
        // Source location
        double x0 = 5.0;
        double y0 = 5.0;
        
        // Create isotropic velocity model
        double velocity = 2.0;  // constant velocity
        
        // Create coefficient matrices for isotropic case
        Eigen::MatrixXd a = Eigen::MatrixXd::Constant(nx, ny, 1.0/(velocity*velocity));
        Eigen::MatrixXd b = Eigen::MatrixXd::Constant(nx, ny, 1.0/(velocity*velocity));
        Eigen::MatrixXd c = Eigen::MatrixXd::Zero(nx, ny);  // no anisotropy
        Eigen::MatrixXd fun = Eigen::MatrixXd::Ones(nx, ny);
        
        // Test CPU computation
        std::cout << "\n=== Running CPU solver ===" << std::endl;
        Eigen::MatrixXd T_cpu;
        auto start_cpu = std::chrono::high_resolution_clock::now();
        EikonSolver2D::FSM_WENO3_PS_2d(xx, yy, a, b, c, fun, x0, y0, T_cpu, false);
        auto end_cpu = std::chrono::high_resolution_clock::now();
        auto duration_cpu = std::chrono::duration_cast<std::chrono::milliseconds>(end_cpu - start_cpu);
        
#ifdef CUDA_ENABLED
        // Test GPU computation
        std::cout << "\n=== Running GPU solver ===" << std::endl;
        Eigen::MatrixXd T_gpu;
        auto start_gpu = std::chrono::high_resolution_clock::now();
        EikonSolver2D::FSM_WENO3_PS_2d(xx, yy, a, b, c, fun, x0, y0, T_gpu, true);
        auto end_gpu = std::chrono::high_resolution_clock::now();
        auto duration_gpu = std::chrono::duration_cast<std::chrono::milliseconds>(end_gpu - start_gpu);
        
        // Compare CPU and GPU results
        double max_diff_cpu_gpu = (T_cpu - T_gpu).cwiseAbs().maxCoeff();
        double mean_diff_cpu_gpu = (T_cpu - T_gpu).cwiseAbs().mean();
        
        std::cout << "\n=== Performance Comparison ===" << std::endl;
        std::cout << "CPU time: " << duration_cpu.count() << " ms" << std::endl;
        std::cout << "GPU time: " << duration_gpu.count() << " ms" << std::endl;
        std::cout << "Speedup: " << (double)duration_cpu.count() / duration_gpu.count() << "x" << std::endl;
        
        std::cout << "\n=== CPU vs GPU Results ===" << std::endl;
        std::cout << "Max difference: " << max_diff_cpu_gpu << std::endl;
        std::cout << "Mean difference: " << mean_diff_cpu_gpu << std::endl;
        
        if (max_diff_cpu_gpu < 1e-3) {
            std::cout << "✓ CPU and GPU results are consistent" << std::endl;
        } else {
            std::cout << "! CPU and GPU results differ significantly" << std::endl;
        }
        
        // Use GPU result for accuracy check
        Eigen::MatrixXd& T_result = T_gpu;
        std::string mode = "GPU";
#else
        std::cout << "\nCUDA support not enabled. Only CPU computation available." << std::endl;
        Eigen::MatrixXd& T_result = T_cpu;
        std::string mode = "CPU";
#endif
        
        
        // Save results
        HDF5Utils hdf5_utils;
        try {
            std::cout << "\n=== Saving results ===" << std::endl;
#ifdef CUDA_ENABLED
            hdf5_utils.saveTravelTimeResults(xx, yy, T_gpu, "travel_time_gpu_result.h5");
            hdf5_utils.saveTravelTimeResults(xx, yy, T_cpu, "travel_time_cpu_result.h5");
#else
            hdf5_utils.saveTravelTimeResults(xx, yy, T_cpu, "travel_time_cpu_result.h5");
#endif
            std::cout << "Results saved to HDF5 files" << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "HDF5 save error: " << e.what() << std::endl;
        }
        
        // Print sample values
        std::cout << "\n=== Sample Values ===" << std::endl;
        int mid_x = nx / 2;
        int mid_y = ny / 2;
        std::cout << "At grid center (" << xx(mid_x) << ", " << yy(mid_y) << "):" << std::endl;
        std::cout << "  CPU: " << T_cpu(mid_x, mid_y) << std::endl;
#ifdef CUDA_ENABLED
        std::cout << "  GPU: " << T_gpu(mid_x, mid_y) << std::endl;
#endif
        
        std::cout << "\n=== Test completed successfully ===" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    }
    
    return 0;
}
