# Eikonal Solver 2D - C++ Implementation

## Overview

This is a C++ rewrite of the original Fortran function `FSM_WENO3_PS_2d`, using the Eigen library for vectorized operations. The solver implements the Fast Sweeping Method combined with WENO3 (Weighted Essentially Non-Oscillatory) scheme to solve the 2D Eikonal equation.

## Function Description

### Original Fortran Function Analysis

The original `FSM_WENO3_PS_2d` function is a numerical solver for the Eikonal equation in 2D anisotropic media:

```math
|∇T|² = f²
```

Where:

- `T` is the traveltime field
- `f` is the slowness field (reciprocal of velocity)
- The equation in anisotropic media is more complex, involving coefficient matrices a, b, c

### Key Features

1. **WENO3 Scheme**: Uses third-order WENO scheme for derivative calculation, providing high accuracy and stability
2. **Fast Sweeping Method**: Four-directional sweeping ensures causality of the numerical solution
3. **Anisotropic Support**: Supports anisotropic media through coefficient matrices a, b, c
4. **Boundary Condition Handling**: Automatic boundary condition processing

## C++ Implementation Features

### Libraries Used
- **Eigen**: For matrix and vector operations, providing optimized linear algebra computations
- **STL**: Standard C++ library for basic functionality

### Vectorization Optimization

1. **Matrix Operations**: Uses Eigen's matrix types instead of primitive arrays
2. **Batch Computation**: Uses Eigen's vectorized operations where possible
3. **Memory Management**: Uses Eigen's automatic memory management

### Code Structure
```cpp
class EikonSolver2D {
public:
    static void FSM_WENO3_PS_2d(...);  // Main solver function
    
private:
    // Helper functions
    static void computeT0AndDerivatives(...);
    static void initializeTauAndChange(...);
    static void computeWENODerivatives(...);
    static void updateBoundaryConditions(...);
    static void computeConvergenceMetrics(...);
};
```

## Function Parameters

### Input Parameters
- `xx`: x-coordinate vector (nx elements)
- `yy`: y-coordinate vector (ny elements)
- `a, b, c`: Anisotropic coefficient matrices (nx × ny)
- `fun`: Slowness field (nx × ny)
- `x0, y0`: Source coordinates
- `u`: Reference solution (nx × ny) - only used for error calculation

### Output Parameters
- `T`: Traveltime field (nx × ny)

## Compilation and Usage

### Dependencies
```bash
mamba create -n EikonPP
mamba activate EikonPP
mamba install -c conda-forge eigen cmake c++_compiler hdf5
```

### Compilation
```bash
mkdir build
cd build
cmake ..
make
```

### Usage Example
```cpp
#include "eikon_solver_2d.hpp"

// Create grid and parameters
Eigen::VectorXd xx = Eigen::VectorXd::LinSpaced(nx, x_min, x_max);
Eigen::VectorXd yy = Eigen::VectorXd::LinSpaced(ny, y_min, y_max);

// Initialize coefficient matrices
Eigen::MatrixXd a(nx, ny), b(nx, ny), c(nx, ny);
Eigen::MatrixXd fun(nx, ny);  // Slowness field
Eigen::MatrixXd u(nx, ny);    // Reference solution

// Set parameters...

// Solve
Eigen::MatrixXd T(nx, ny);
EikonSolver2D::FSM_WENO3_PS_2d(xx, yy, a, b, c, fun, x0, y0, T);
```

## Algorithm Details

### 1. Initialization Phase
- Calculate grid spacing
- Perform bilinear interpolation at source location to obtain local parameters
- Compute initial traveltime field T0 and its derivatives
- Initialize tau array and change flag array

### 2. Iterative Solution
- Four-directional sweeping (forward and backward in both directions)
- For each point to be updated:
  - Compute WENO3 derivatives
  - Calculate Hamiltonian Htau
  - Update tau values
- Update boundary conditions
- Check convergence

### 3. WENO3 Scheme
Uses third-order Weighted Essentially Non-Oscillatory scheme:
- Calculate weights based on local smoothness indicators
- Combine different finite difference schemes
- Ensure stability and accuracy of the numerical solution

## Performance Optimization

### Advantages of C++ over Fortran
1. **Memory Management**: Eigen provides efficient memory management
2. **Vectorization**: Modern compilers can better vectorize C++ code
3. **Template Optimization**: Compile-time optimization
4. **SIMD Support**: Eigen automatically uses SIMD instructions

### Performance Tuning Recommendations
1. Use `-O3 -march=native` compilation flags
2. Consider OpenMP parallelization (can be added in loops)
3. Use appropriate data types (double vs float)

## Differences from Original Fortran Version

### Main Improvements
1. **Type Safety**: C++'s strong type system
2. **Memory Safety**: Automatic memory management
3. **Interface Simplification**: Clear class interface
4. **Extensibility**: Easy to add new features

### Numerical Precision
- Maintains the same numerical precision as the original Fortran version
- Same convergence criteria and algorithm logic

## Testing and Validation

Recommended test cases:
1. Analytical solution verification in uniform media
2. Anisotropic media testing
3. Complex topography testing
4. Performance benchmarking

## Extension Possibilities

1. **Parallelization**: OpenMP or CUDA implementation
2. **3D Extension**: Extend to three-dimensional problems
3. **Adaptive Mesh**: Implement adaptive mesh refinement
4. **Other Schemes**: Implement other high-order schemes (WENO5, etc.)

## Troubleshooting

Common issues:
1. **Compilation Errors**: Check Eigen installation
2. **Numerical Instability**: Check grid resolution and parameter settings
3. **Convergence Problems**: Adjust tolerance parameters

This C++ implementation maintains the numerical precision and algorithm logic of the original Fortran code while providing the advantages of modern C++ and the performance optimizations of the Eigen library.
