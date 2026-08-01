# EiknoalPrism

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++](https://img.shields.io/badge/C++-17-blue.svg)](https://isocpp.org/)
[![CMake](https://img.shields.io/badge/CMake-3.10+-green.svg)](https://cmake.org/)

A modern C++ implementation for processing Source Functions (SF) using the Eikonal equation to generate distance fields and prism meshes.

## 📋 Table of Contents

- [Features](#-features)
- [Installation](#-installation)
- [Pipeline (What It Does)](#-Pipeline)
- [Usage](#-usage)
- [Output Files](#-output-files)
- [Dependencies](#-dependencies)
- [License](#-license)

## ✨ Features

- **Source Function Processing**: Read closed polygon representations of Source Functions
- **Background Cartesian Grid**: Automatic grid generation based on input geometry
- **Eikonal Equation Solver**: Compute accurate distance fields using the Eikonal equation
- **Gradient Computation**: Calculate gradient magnitude of distance functions
- **Prism Mesh Generation**: Create high-quality prism meshes for Source Functions
- **Visualization Ready**: Export VTK and MSH formats for visualization

## 🔧 Installation

### Prerequisites

- C++17 compatible compiler (GCC 7+, Clang 5+, or MSVC 2017+)
- CMake 3.10 or higher

### Build Instructions

```bash
# Clone the repository
git clone git@github.com:graceidea/eiknoalPrism.git
cd eiknoalPrism/

# Build the binaries
mkdir build && cd build
cmake ..
make -j$(nproc)  # Use all available cores
ls ../../tests/example/

# Expected output:
  input.obj polygon.msh grid.msh distance_exact.vtk gradient_magnitude.vtk prism.msh

### 🧪 Running Tests
cd eiknoalPrism/build/bin/
./test


### ⚙️ Pipeline (What It Does)

The `eikonalPrism` engine executes a 5-step geometric and differential pipeline to generate prism meshes from a closed polygon boundary using the Eikonal equation.

| Step | Operation Phase | Input Asset | Output Generated File | Output Visualization |
| :---: | :--- | :--- | :--- | :---: |
| **1** | **Source Function (SF) Processing** | `input.obj` | `polygon.msh` | [View Geometry](#1-source-function) |
| **2** | **Background Cartesian Grid (BCG)** | `polygon.msh` | `grid.msh` | [View Grid](#2-background-cartesian-grid) |
| **3** | **Eikonal Distance Field Computation** | `grid.msh` | `distance_exact.vtk` | [View Field](#3-distance-field) |
| **4** | **Gradient of Distance Function (GDF)**| `distance_exact.vtk` | `gradient_magnitude.vtk`| [View Gradient](#4-gradient-field) |
| **5** | **Prism Mesh Generation** | `gradient_magnitude.vtk`| `prism.msh` | [View Mesh](#5-prism-mesh) |

---



The application processes Source Functions (SF) represented as closed polygons through a five-step pipeline:

1. Read Source Function (Input)
The process begins by reading a closed polygon representation of the Source Function from an OBJ file. The input file is located at eiknoalPrism/tests/example/input.obj.

<div align="center"> <img width="573" height="567" alt="Input Source Function" src="https://github.com/user-attachments/assets/60c52924-c825-404e-9e0e-464f4000b77d" /> <br> <em><strong>Figure 1:</strong> Input Source Function (SF) - Closed polygon representation</em> </div>
2. Generate Background Cartesian Grid (BCG)
A structured Cartesian grid is automatically generated based on the input Source Function boundaries. This grid provides the computational domain for subsequent calculations. The grid file is saved as grid.msh in the tests/example/ directory.

<div align="center"> <img width="547" height="550" alt="BCG and SF" src="https://github.com/user-attachments/assets/3de02070-47e9-49e8-b4d7-212434b02843" /> <br> <em><strong>Figure 2:</strong> Background Cartesian Grid (BCG) with Source Function overlay</em> </div>
3. Compute Distance Function
The Eikonal equation is solved to compute the signed distance field from the Source Function. This distance field is crucial for understanding the geometric properties of the Source Function. The output is saved as distance_exact.vtk.

<div align="center"> <img width="530" height="532" alt="Distance Field" src="https://github.com/user-attachments/assets/ba4d7b17-d050-4af9-8cf4-bfb7cf9b7eb0" /> <br> <em><strong>Figure 3:</strong> Distance field computed from the Source Function using the Eikonal equation</em> </div>
4. Compute Gradient Magnitude
The gradient magnitude of the distance function is calculated to analyze geometric properties such as boundary normals and curvature. This output is saved as gradient_magnitude.vtk.

<div align="center"> <img width="529" height="533" alt="Gradient Magnitude" src="https://github.com/user-attachments/assets/bbf6c3c2-ff62-4d83-ac45-c895d6ec39a6" /> <br> <em><strong>Figure 4:</strong> Gradient magnitude of the distance function</em> </div>
5. Generate Prism Mesh
The final step generates a high-quality prism mesh for the Source Function with proper boundary layer resolution. This mesh is ideal for finite element simulations and is saved as prism.msh.

<div align="center"> <img width="505" height="534" alt="Prism Mesh" src="https://github.com/user-attachments/assets/2c1f973c-178f-4e40-9a9a-59cc7c8848ab" /> <br> <em><strong>Figure 5:</strong> Prism mesh generated for the Source Function</em> </div>
📂 Usage
Running the Test Suite
bash
cd build/bin/
./test
Input/Output Files
Stage	File	Location	Format	Description
Input	input.obj	tests/example/	OBJ	Source Function as closed polygon
Output	polygon.msh	tests/example/	MSH	Mesh representation of input polygon
Output	grid.msh	tests/example/	MSH	Background Cartesian grid
Output	distance_exact.vtk	tests/example/	VTK	Exact distance field (Eikonal solution)
Output	gradient_magnitude.vtk	tests/example/	VTK	Gradient magnitude of distance field
Output	prism.msh	tests/example/	MSH	Final prism mesh
Example Verification
After running ./test, you can verify all output files are generated:

bash
cd build/bin/
./test


📊 Visualization
ParaView (Recommended for .vtk files) - Open-source scientific visualization
Gmsh (Recommended for polygon.msh, grid.msh, prism.msh) - Mesh generation and visualization



🔬 Dependencies
HDF5	Mesh file I/O	1.10+
Eigen3	Linear algebra operations	3.3+
Boost	Utility functions	1.70+

📄 License
This project is licensed under the MIT License - see the LICENSE file for details.












