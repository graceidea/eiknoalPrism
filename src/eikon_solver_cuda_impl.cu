#ifdef CUDA_ENABLED
#include "eikon_solver_2d.hpp"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <thrust/extrema.h>
#include <thrust/transform_reduce.h>
#include <thrust/device_ptr.h>
#include <thrust/functional.h>
#include <thrust/reduce.h>
#include <stdexcept>
#include <cmath>

// Static GPU memory pointers initialization
double* EikonSolver2D::d_tau = nullptr;
double* EikonSolver2D::d_tau_old = nullptr;
double* EikonSolver2D::d_T0 = nullptr;
double* EikonSolver2D::d_T0x = nullptr;
double* EikonSolver2D::d_T0y = nullptr;
double* EikonSolver2D::d_a = nullptr;
double* EikonSolver2D::d_b = nullptr;
double* EikonSolver2D::d_c = nullptr;
double* EikonSolver2D::d_fun = nullptr;
int* EikonSolver2D::d_ischange = nullptr;
double* EikonSolver2D::d_sigx_matrix = nullptr;
double* EikonSolver2D::d_sigy_matrix = nullptr;
double* EikonSolver2D::d_inv_dx_dy = nullptr;
double* EikonSolver2D::d_H_const = nullptr;
double* EikonSolver2D::d_H_coeff_x = nullptr;
double* EikonSolver2D::d_H_coeff_y = nullptr;
double* EikonSolver2D::d_reduction_buffer = nullptr;
double* EikonSolver2D::d_diff_buffer = nullptr; 
double* EikonSolver2D::d_T = nullptr;
bool EikonSolver2D::gpu_initialized = false;

__device__ const int blockx = 32;
__device__ const int blocky = 16;

// Device math functions using CUDA built-in functions
__device__  double device_pow(double x, double y) {
    return pow(x, y);  // Use CUDA intrinsic
}

__device__  double device_fabs(double x) {
    return fabs(x);    // Use CUDA intrinsic
}

__device__  double device_fmax(double x, double y) {
    return fmax(x, y); // Use CUDA intrinsic
}

__device__  double device_sqrt(double x) {
    return sqrt(x); // Use CUDA intrinsic with round-to-nearest
}

// Helper function to convert 2D indices to 1D index for Eigen column-major matrices
// Eigen matrices are stored in column-major order: M(i,j) = data[j*rows + i]
// For our case: tau(ix, iy) = data[iy*nx + ix]
__device__ inline int I2V_cuda(int ix, int iy, int nx) {
    return iy * nx + ix;
}

// WENO3 derivatives computation function
__device__ void computeWENO3Derivatives(
    const double* __restrict__ tau,
    int iix, int iiy, int nx, int ny,
    double dx, double dy,
    double& px1, double& px2, double& py1, double& py2
) {
    const double eps = 1e-14;
    px1 = px2 = py1 = py2 = 0.0;
    
    // X-direction derivatives
    if (iix == 1) {
        px1 = (tau[I2V_cuda(iix, iiy, nx)] - tau[I2V_cuda(iix-1, iiy, nx)]) / dx;
        
        double numerator = eps + device_pow(tau[I2V_cuda(iix, iiy, nx)] - 2 * tau[I2V_cuda(iix+1, iiy, nx)] + tau[I2V_cuda(iix+2, iiy, nx)], 2.0);
        double denominator = eps + device_pow(tau[I2V_cuda(iix-1, iiy, nx)] - 2 * tau[I2V_cuda(iix, iiy, nx)] + tau[I2V_cuda(iix+1, iiy, nx)], 2.0);
        double wx2 = 1.0 / (1.0 + 2.0 * device_pow(numerator / denominator, 2.0));
        
        px2 = (1 - wx2) * (tau[I2V_cuda(iix+1, iiy, nx)] - tau[I2V_cuda(iix-1, iiy, nx)]) / (2.0 * dx) +
              wx2 * (-3 * tau[I2V_cuda(iix, iiy, nx)] + 4 * tau[I2V_cuda(iix+1, iiy, nx)] - tau[I2V_cuda(iix+2, iiy, nx)]) / (2.0 * dx);
    }
    else if (iix == nx - 2) {
        double numerator = eps + device_pow(tau[I2V_cuda(iix, iiy, nx)] - 2 * tau[I2V_cuda(iix-1, iiy, nx)] + tau[I2V_cuda(iix-2, iiy, nx)], 2.0);
        double denominator = eps + device_pow(tau[I2V_cuda(iix+1, iiy, nx)] - 2 * tau[I2V_cuda(iix, iiy, nx)] + tau[I2V_cuda(iix-1, iiy, nx)], 2.0);
        double wx1 = 1.0 / (1.0 + 2.0 * device_pow(numerator / denominator, 2.0));
        
        px1 = (1 - wx1) * (tau[I2V_cuda(iix+1, iiy, nx)] - tau[I2V_cuda(iix-1, iiy, nx)]) / (2.0 * dx) +
              wx1 * (3 * tau[I2V_cuda(iix, iiy, nx)] - 4 * tau[I2V_cuda(iix-1, iiy, nx)] + tau[I2V_cuda(iix-2, iiy, nx)]) / (2.0 * dx);
        px2 = (tau[I2V_cuda(iix+1, iiy, nx)] - tau[I2V_cuda(iix, iiy, nx)]) / dx;
    }
    else if (iix >= 2 && iix < nx - 2) {
        // Interior points
        double numerator1 = eps + device_pow(tau[I2V_cuda(iix, iiy, nx)] - 2 * tau[I2V_cuda(iix-1, iiy, nx)] + tau[I2V_cuda(iix-2, iiy, nx)], 2.0);
        double denominator1 = eps + device_pow(tau[I2V_cuda(iix+1, iiy, nx)] - 2 * tau[I2V_cuda(iix, iiy, nx)] + tau[I2V_cuda(iix-1, iiy, nx)], 2.0);
        double wx1 = 1.0 / (1.0 + 2.0 * device_pow(numerator1 / denominator1, 2.0));

        px1 = (1.0 - wx1) * (tau[I2V_cuda(iix+1, iiy, nx)] - tau[I2V_cuda(iix-1, iiy, nx)]) / (2.0 * dx) +
              wx1 * (3 * tau[I2V_cuda(iix, iiy, nx)] - 4 * tau[I2V_cuda(iix-1, iiy, nx)] + tau[I2V_cuda(iix-2, iiy, nx)]) / (2.0 * dx);

        double numerator2 = eps + device_pow(tau[I2V_cuda(iix, iiy, nx)] - 2 * tau[I2V_cuda(iix+1, iiy, nx)] + tau[I2V_cuda(iix+2, iiy, nx)], 2.0);
        double denominator2 = eps + device_pow(tau[I2V_cuda(iix-1, iiy, nx)] - 2 * tau[I2V_cuda(iix, iiy, nx)] + tau[I2V_cuda(iix+1, iiy, nx)], 2.0);
        double wx2 = 1.0 / (1.0 + 2.0 * device_pow(numerator2 / denominator2, 2.0));

        px2 = (1.0 - wx2) * (tau[I2V_cuda(iix+1, iiy, nx)] - tau[I2V_cuda(iix-1, iiy, nx)]) / (2.0 * dx) +
              wx2 * (-3 * tau[I2V_cuda(iix, iiy, nx)] + 4 * tau[I2V_cuda(iix+1, iiy, nx)] - tau[I2V_cuda(iix+2, iiy, nx)]) / (2.0 * dx);
    }
    
    // Y-direction derivatives  
    if (iiy == 1) {
        py1 = (tau[I2V_cuda(iix, iiy, nx)] - tau[I2V_cuda(iix, iiy-1, nx)]) / dy;

        double numerator = eps + device_pow(tau[I2V_cuda(iix, iiy, nx)] - 2 * tau[I2V_cuda(iix, iiy+1, nx)] + tau[I2V_cuda(iix, iiy+2, nx)], 2.0);
        double denominator = eps + device_pow(tau[I2V_cuda(iix, iiy-1, nx)] - 2 * tau[I2V_cuda(iix, iiy, nx)] + tau[I2V_cuda(iix, iiy+1, nx)], 2.0);
        double wy2 = 1.0 / (1.0 + 2.0 * device_pow(numerator / denominator, 2.0));

        py2 = (1 - wy2) * (tau[I2V_cuda(iix, iiy+1, nx)] - tau[I2V_cuda(iix, iiy-1, nx)]) / (2.0 * dy) +
              wy2 * (-3 * tau[I2V_cuda(iix, iiy, nx)] + 4 * tau[I2V_cuda(iix, iiy+1, nx)] - tau[I2V_cuda(iix, iiy+2, nx)]) / (2.0 * dy);
    }
    else if (iiy == ny - 2) {
        double numerator = eps + device_pow(tau[I2V_cuda(iix, iiy, nx)] - 2 * tau[I2V_cuda(iix, iiy-1, nx)] + tau[I2V_cuda(iix, iiy-2, nx)], 2.0);
        double denominator = eps + device_pow(tau[I2V_cuda(iix, iiy+1, nx)] - 2 * tau[I2V_cuda(iix, iiy, nx)] + tau[I2V_cuda(iix, iiy-1, nx)], 2.0);
        double wy1 = 1.0 / (1.0 + 2.0 * device_pow(numerator / denominator, 2.0));

        py1 = (1 - wy1) * (tau[I2V_cuda(iix, iiy+1, nx)] - tau[I2V_cuda(iix, iiy-1, nx)]) / (2.0 * dy) +
              wy1 * (3 * tau[I2V_cuda(iix, iiy, nx)] - 4 * tau[I2V_cuda(iix, iiy-1, nx)] + tau[I2V_cuda(iix, iiy-2, nx)]) / (2.0 * dy);
        py2 = (tau[I2V_cuda(iix, iiy+1, nx)] - tau[I2V_cuda(iix, iiy, nx)]) / dy;
    }
    else if (iiy >= 2 && iiy < ny - 2) {
        // Interior points
        double numerator1 = eps + device_pow(tau[I2V_cuda(iix, iiy, nx)] - 2 * tau[I2V_cuda(iix, iiy-1, nx)] + tau[I2V_cuda(iix, iiy-2, nx)], 2.0);
        double denominator1 = eps + device_pow(tau[I2V_cuda(iix, iiy+1, nx)] - 2 * tau[I2V_cuda(iix, iiy, nx)] + tau[I2V_cuda(iix, iiy-1, nx)], 2.0);
        double wy1 = 1.0 / (1.0 + 2.0 * device_pow(numerator1 / denominator1, 2.0));

        py1 = (1 - wy1) * (tau[I2V_cuda(iix, iiy+1, nx)] - tau[I2V_cuda(iix, iiy-1, nx)]) / (2.0 * dy) +
              wy1 * (3 * tau[I2V_cuda(iix, iiy, nx)] - 4 * tau[I2V_cuda(iix, iiy-1, nx)] + tau[I2V_cuda(iix, iiy-2, nx)]) / (2.0 * dy);

        double numerator2 = eps + device_pow(tau[I2V_cuda(iix, iiy, nx)] - 2 * tau[I2V_cuda(iix, iiy+1, nx)] + tau[I2V_cuda(iix, iiy+2, nx)], 2.0);
        double denominator2 = eps + device_pow(tau[I2V_cuda(iix, iiy-1, nx)] - 2 * tau[I2V_cuda(iix, iiy, nx)] + tau[I2V_cuda(iix, iiy+1, nx)], 2.0);
        double wy2 = 1.0 / (1.0 + 2.0 * device_pow(numerator2 / denominator2, 2.0));

        py2 = (1 - wy2) * (tau[I2V_cuda(iix, iiy+1, nx)] - tau[I2V_cuda(iix, iiy-1, nx)]) / (2.0 * dy) +
              wy2 * (-3 * tau[I2V_cuda(iix, iiy, nx)] + 4 * tau[I2V_cuda(iix, iiy+1, nx)] - tau[I2V_cuda(iix, iiy+2, nx)]) / (2.0 * dy);
    }
}

// 简化的差值计算kernel - 计算tau和tau_old的差值，同时支持L1和Linf范数
__global__ void computeDifferenceKernel(
    const double* __restrict__ tau,
    const double* __restrict__ tau_old,
    double* __restrict__ diff_buffer,
    int nx, int ny,
    double dx, double dy
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (idx < nx && idy < ny) {
        int id = I2V_cuda(idx, idy, nx);
        double diff = device_fabs(tau[id] - tau_old[id]);
        // 存储差值的绝对值，用于后续的L1和Linf计算
        diff_buffer[id] = diff;
    }
}

// Copy tau to tau_old on GPU
__global__ void copyTauToOldKernel(const double* __restrict__ tau, double* __restrict__ tau_old, const int nx, const int ny) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;
    int id = I2V_cuda(idx, idy, nx);
    if (id < nx * ny) {
        tau_old[id] = tau[id];
    }
}

// GPU boundary conditions update kernel
__global__ void updateBoundaryConditionsGPU(double* __restrict__ tau, int nx, int ny) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;
    
    // Left and right boundaries
    if (idx == 0 && idy < ny) {
        tau[I2V_cuda(0, idy, nx)] = device_fmax(2.0 * tau[I2V_cuda(1, idy, nx)] - tau[I2V_cuda(2, idy, nx)],
                                         tau[I2V_cuda(2, idy, nx)]);
    }
    if (idx == nx-1 && idy < ny) {
        tau[I2V_cuda(nx-1, idy, nx)] = device_fmax(2.0 * tau[I2V_cuda(nx-2, idy, nx)] - tau[I2V_cuda(nx-3, idy, nx)],
                                            tau[I2V_cuda(nx-3, idy, nx)]);
    }
    
    // Top and bottom boundaries
    if (idy == 0 && idx < nx) {
        tau[I2V_cuda(idx, 0, nx)] = device_fmax(2.0 * tau[I2V_cuda(idx, 1, nx)] - tau[I2V_cuda(idx, 2, nx)],
                                         tau[I2V_cuda(idx, 2, nx)]);
    }
    if (idy == ny-1 && idx < nx) {
        tau[I2V_cuda(idx, ny-1, nx)] = device_fmax(2.0 * tau[I2V_cuda(idx, ny-2, nx)] - tau[I2V_cuda(idx, ny-3, nx)],
                                            tau[I2V_cuda(idx, ny-3, nx)]);
    }
}

// CUDA kernel for GPU sweep computation
__global__ void gpu_sweep_kernel(
    double* __restrict__ tau,
    const double* __restrict__ a,
    const double* __restrict__ b,
    const double* __restrict__ c,
    const double* __restrict__ fun,
    const double* __restrict__ T0x,
    const double* __restrict__ T0y,
    const int* __restrict__ ischange,
    const double* __restrict__ sigx_matrix,
    const double* __restrict__ sigy_matrix,
    const double* __restrict__ inv_dx_dy,
    const double* __restrict__ H_const,
    const double* __restrict__ H_coeff_x,
    const double* __restrict__ H_coeff_y,
    const int x_start, const int x_end, 
    const int y_start, const int y_end,
    const int xdirec, const int ydirec, 
    const double dx, const double dy, 
    const int nx, const int ny
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;
    
    // Map thread indices to sweep coordinates
    // Check bounds
    if (xdirec < 0 ){
        if (idx < x_start || idx >= x_end) return;
    } else {
        if (idx > x_start || idx <= x_end) return;
    }

    if (ydirec < 0) {
        if (idy < y_start || idy >= y_end) return;
    } else {
        if (idy > y_start || idy <= y_end) return;
    }

    // Use I2V function for Eigen column-major indexing
    int id = I2V_cuda(idx, idy, nx);

    if (ischange[id] == 1) {
        // Use pre-computed values
        double sigx = sigx_matrix[id];
        double sigy = sigy_matrix[id];
        double coe = inv_dx_dy[id];
        
        // Compute WENO3 derivatives using separate function
        double px1, px2, py1, py2;
        computeWENO3Derivatives(tau, idx, idy, nx, ny, dx, dy, px1, px2, py1, py2);

        // Optimized Hamiltonian computation
        double px_avg = (px1 + px2) * 0.5;
        double py_avg = (py1 + py2) * 0.5;
        
        double H_quadratic = a[id] * px_avg * px_avg + 
                             b[id] * py_avg * py_avg - 
                             2.0 * c[id] * px_avg * py_avg;
        
        double H_linear = 2.0 * (H_coeff_x[id] * px_avg + H_coeff_y[id] * py_avg);
        
        double Htau = device_sqrt(H_quadratic + H_linear + H_const[id]);
        
        // Update tau (using column-major indexing)
        double derivative_correction = coe * (sigx * (px2 - px1) + sigy * (py2 - py1)) * 0.5;
        double tpT = coe * (fun[id] - Htau) + derivative_correction + tau[id];
        
        tau[id] = tpT;
    }
}

// Initialize GPU memory
void EikonSolver2D::initializeGPU(int nx, int ny) {
    if (gpu_initialized) return;
    
    size_t size = nx * ny * sizeof(double);
    size_t size_int = nx * ny * sizeof(int);
    
    // Check CUDA device
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount == 0) {
        throw std::runtime_error("No CUDA-capable devices found");
    }
    
    // Allocate device memory
    cudaMalloc(&d_tau, size);
    cudaMalloc(&d_tau_old, size);
    cudaMalloc(&d_T0, size);
    cudaMalloc(&d_T0x, size);
    cudaMalloc(&d_T0y, size);
    cudaMalloc(&d_a, size);
    cudaMalloc(&d_b, size);
    cudaMalloc(&d_c, size);
    cudaMalloc(&d_fun, size);
    cudaMalloc(&d_ischange, size_int);
    cudaMalloc(&d_sigx_matrix, size);
    cudaMalloc(&d_sigy_matrix, size);
    cudaMalloc(&d_inv_dx_dy, size);
    cudaMalloc(&d_H_const, size);
    cudaMalloc(&d_H_coeff_x, size);
    cudaMalloc(&d_H_coeff_y, size);
    cudaMalloc(&d_T, size);
    cudaMalloc(&d_diff_buffer, size); 

    cudaMalloc(&d_reduction_buffer, sizeof(double));
    
    gpu_initialized = true;
}

// Copy data to GPU
void EikonSolver2D::copyDataToGPU(
    const Eigen::MatrixXd& a, const Eigen::MatrixXd& b, 
    const Eigen::MatrixXd& c, const Eigen::MatrixXd& fun,
    const Eigen::MatrixXd& T0, const Eigen::MatrixXd& T0x,
    const Eigen::MatrixXd& T0y, const Eigen::MatrixXi& ischange,
    const Eigen::MatrixXd& tau,
    double dx, double dy
) {
    int nx = a.rows();
    int ny = a.cols();
    size_t size = nx * ny * sizeof(double);
    size_t size_int = nx * ny * sizeof(int);
    
    // Copy input matrices
    cudaMemcpy(d_a, a.data(), size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b.data(), size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_c, c.data(), size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_fun, fun.data(), size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_T0, T0.data(), size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_T0x, T0x.data(), size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_T0y, T0y.data(), size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_ischange, ischange.data(), size_int, cudaMemcpyHostToDevice);
    cudaMemcpy(d_tau, tau.data(), size, cudaMemcpyHostToDevice);

    // Pre-compute coefficient matrices on CPU (more reliable than GPU kernel)
    Eigen::MatrixXd sigx_matrix = a.array().sqrt();
    Eigen::MatrixXd sigy_matrix = b.array().sqrt();
    Eigen::MatrixXd inv_dx_dy = (sigx_matrix / dx + sigy_matrix / dy).cwiseInverse();
    
    Eigen::MatrixXd H_const = a.cwiseProduct(T0x.cwiseProduct(T0x)) + 
                              b.cwiseProduct(T0y.cwiseProduct(T0y)) - 
                              2.0 * c.cwiseProduct(T0x.cwiseProduct(T0y));
    
    Eigen::MatrixXd H_coeff_x = a.cwiseProduct(T0x) - c.cwiseProduct(T0y);
    Eigen::MatrixXd H_coeff_y = b.cwiseProduct(T0y) - c.cwiseProduct(T0x);
    
    // Copy pre-computed matrices
    cudaMemcpy(d_sigx_matrix, sigx_matrix.data(), size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_sigy_matrix, sigy_matrix.data(), size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_inv_dx_dy, inv_dx_dy.data(), size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_H_const, H_const.data(), size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_H_coeff_x, H_coeff_x.data(), size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_H_coeff_y, H_coeff_y.data(), size, cudaMemcpyHostToDevice);
}

__global__ void addTauToT0GPU(double* d_tau, double* d_T0, double* d_T, int nx, int ny) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (idx < nx && idy < ny) {
        int id = I2V_cuda(idx, idy, nx);
        d_T[id] = d_tau[id] + d_T0[id];
    }
}

// Optimized batch sweep to reduce memory transfer overhead
void EikonSolver2D::sweepAllDirectionsGPU(
    Eigen::MatrixXd& T, double dx, double dy, int nx, int ny
) {
    double L1_dif, Linf_dif;

    // Initialize tau on GPU (copy from T0)
    int nlen = nx * ny;
    size_t size = nlen * sizeof(double);
    
    dim3 blockSize(blockx, blocky);
    int grid_x = (nx + blockSize.x - 1) / blockSize.x;
    int grid_y = (ny + blockSize.y - 1) / blockSize.y;
    dim3 gridSize(grid_x, grid_y);

    for (int iter = 0; iter < MaxIter; ++iter) {
        copyTauToOldKernel<<<gridSize, blockSize>>>(d_tau, d_tau_old, nx, ny);

        for (int xdirec = -1; xdirec <= 1; xdirec += 2) {
            for (int ydirec = -1; ydirec <= 1; ydirec += 2) {

                int x_start, x_end, y_start, y_end;
                selectDirection(xdirec, ydirec, nx, ny, x_start, x_end, y_start, y_end);
                
                gpu_sweep_kernel<<<gridSize, blockSize>>>(
                    d_tau, d_a, d_b, d_c, d_fun, d_T0x, d_T0y, d_ischange,
                    d_sigx_matrix, d_sigy_matrix, d_inv_dx_dy,
                    d_H_const, d_H_coeff_x, d_H_coeff_y,
                    x_start, x_end, y_start, y_end,
                    xdirec, ydirec, dx, dy, nx, ny
                );
                
                cudaDeviceSynchronize();
                
                // Update boundary conditions on GPU after each sweep (like CPU version)
                updateBoundaryConditionsGPU<<<gridSize, blockSize>>>(d_tau, nx, ny);
                cudaDeviceSynchronize();
            }
        }
        computeSimpleConvergenceMetricsGPU(nx, ny, dx, dy, L1_dif, Linf_dif);
        if (std::abs(L1_dif) < tol && Linf_dif < tol) {
            break;
        }
    }
    // After all sweeps, compute final T
    addTauToT0GPU<<<gridSize, blockSize>>>(d_tau, d_T0, d_T, nx, ny);
    cudaDeviceSynchronize();

    cudaMemcpy(T.data(), d_T, size, cudaMemcpyDeviceToHost);
    cleanupGPU();
}

// 简化的GPU收敛性检查 - 计算L1_dif和Linf_dif，使用Thrust库
void EikonSolver2D::computeSimpleConvergenceMetricsGPU(
    int nx, int ny, double dx, double dy, double& L1_dif, double& Linf_dif
) {
    // 计算差值
    dim3 blockSize(blockx, blocky);
    dim3 gridSize((nx + blockSize.x - 1) / blockSize.x, (ny + blockSize.y - 1) / blockSize.y);
    
    computeDifferenceKernel<<<gridSize, blockSize>>>(d_tau, d_tau_old, d_diff_buffer, nx, ny, dx, dy);
    cudaDeviceSynchronize();
    
    // 使用Thrust库进行高效的reduction
    thrust::device_ptr<double> diff_ptr(d_diff_buffer);
    
    // 计算Linf范数（最大值）
    Linf_dif = *thrust::max_element(diff_ptr, diff_ptr + nx * ny);
    
    // 计算L1范数（所有差值的和乘以面积元素）
    double sum_diff = thrust::reduce(diff_ptr, diff_ptr + nx * ny, 0.0, thrust::plus<double>());
    L1_dif = sum_diff * dx * dy;
}

// Cleanup GPU memory
void EikonSolver2D::cleanupGPU() {
    if (gpu_initialized) {
        cudaFree(d_tau);
        cudaFree(d_tau_old);
        cudaFree(d_T0);
        cudaFree(d_T0x);
        cudaFree(d_T0y);
        cudaFree(d_a);
        cudaFree(d_b);
        cudaFree(d_c);
        cudaFree(d_fun);
        cudaFree(d_ischange);
        cudaFree(d_sigx_matrix);
        cudaFree(d_sigy_matrix);
        cudaFree(d_inv_dx_dy);
        cudaFree(d_H_const);
        cudaFree(d_H_coeff_x);
        cudaFree(d_H_coeff_y);
        cudaFree(d_reduction_buffer);
        cudaFree(d_diff_buffer);
        cudaFree(d_T);
        
        gpu_initialized = false;
    }
}

#endif // CUDA_ENABLED
