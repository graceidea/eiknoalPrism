#include <iostream>
#include <fstream>      // ← ADD THIS for std::ofstream
#include <string>       // ← ADD THIS for std::string
#include <chrono>
#include <Eigen/Dense>

#include "eikon_solver_2d.hpp"
#include "hdf5_utils.hpp"

#include <cmath>


// Simple structure for triangle
struct Triangle2D {
    double x1, y1, x2, y2, x3, y3;
    
    Triangle2D(double x1_, double y1_, double x2_, double y2_, double x3_, double y3_)
        : x1(x1_), y1(y1_), x2(x2_), y2(y2_), x3(x3_), y3(y3_) {}
};

// Function to compute signed distance to triangle (for reference/validation)
double distanceToTriangle(double x, double y, const Triangle2D& tri) {
    // Simple distance to centroid (for testing)
    double cx = (tri.x1 + tri.x2 + tri.x3) / 3.0;
    double cy = (tri.y1 + tri.y2 + tri.y3) / 3.0;
    return std::sqrt((x - cx)*(x - cx) + (y - cy)*(y - cy));
}

// Function to save as VTK
void saveAsVTK(const Eigen::MatrixXd& T, const Eigen::VectorXd& xx, 
               const Eigen::VectorXd& yy, const std::string& filename,
               const std::string& fieldname) {
    int nx = xx.size();
    int ny = yy.size();
    
    std::ofstream vtk_file(filename);
    if (!vtk_file.is_open()) {
        std::cerr << "Error: Could not open " << filename << std::endl;
        return;
    }
    
    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "Distance Field\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET STRUCTURED_GRID\n";
    vtk_file << "DIMENSIONS " << nx << " " << ny << " 1\n";
    vtk_file << "POINTS " << nx * ny << " float\n";
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            vtk_file << xx(i) << " " << yy(j) << " 0.0\n";
        }
    }
    
    vtk_file << "POINT_DATA " << nx * ny << "\n";
    vtk_file << "SCALARS " << fieldname << " float 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            vtk_file << T(i, j) << "\n";
        }
    }
    
    vtk_file.close();
    std::cout << "VTK file saved: " << filename << std::endl;
}

// Function to save as CSV
void saveAsCSV(const Eigen::MatrixXd& T, const std::string& filename) {
    std::ofstream csv_file(filename);
    if (!csv_file.is_open()) {
        std::cerr << "Error: Could not open " << filename << std::endl;
        return;
    }
    
    for (int i = 0; i < T.rows(); ++i) {
        for (int j = 0; j < T.cols(); ++j) {
            csv_file << T(i, j);
            if (j < T.cols() - 1) csv_file << ",";
        }
        csv_file << "\n";
    }
    
    csv_file.close();
    std::cout << "CSV file saved: " << filename << std::endl;
}
// ========================================================================
// WRITE TRIANGLE AS GMSH FILE
// ========================================================================
void writeTriangleGmsh(const Triangle2D& tri, const std::string& filename) {
    std::ofstream msh_file(filename);
    if (!msh_file.is_open()) {
        std::cerr << "Error: Could not open " << filename << std::endl;
        return;
    }
    
    // Gmsh ASCII format header
    msh_file << "$MeshFormat\n";
    msh_file << "2.2 0 8\n";
    msh_file << "$EndMeshFormat\n";
    
    // Write nodes (3 vertices)
    msh_file << "$Nodes\n";
    msh_file << "3\n";  // 3 nodes
    msh_file << "1 " << tri.x1 << " " << tri.y1 << " 0.0\n";
    msh_file << "2 " << tri.x2 << " " << tri.y2 << " 0.0\n";
    msh_file << "3 " << tri.x3 << " " << tri.y3 << " 0.0\n";
    msh_file << "$EndNodes\n";
    
    // Write elements (1 triangle)
    msh_file << "$Elements\n";
    msh_file << "1\n";  // 1 element
    // Element type 2 = triangle, with 2 tags, connected to nodes 1, 2, 3
    msh_file << "1 2 2 0 0 1 2 3\n";
    msh_file << "$EndElements\n";
    
    // Write node data (optional - color the triangle vertices)
    msh_file << "$NodeData\n";
    msh_file << "1\n";  // 1 field
    msh_file << "\"Triangle\"\n";
    msh_file << "1\n";  // Real data
    msh_file << "0.0\n";  // Time
    msh_file << "3\n";  // 3 components
    msh_file << "3\n";  // 3 data entries
    msh_file << "1 1.0\n";
    msh_file << "2 1.0\n";
    msh_file << "3 1.0\n";
    msh_file << "$EndNodeData\n";
    
    msh_file.close();
    std::cout << "Triangle Gmsh file saved: " << filename << std::endl;
}

int main() {
    std::cout << "==========================================" << std::endl;
    std::cout << "  TRIANGLE DISTANCE FIELD USING EIKONAL SOLVER" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    // Define grid
    int nx = 201, ny = 201;
    double x_min = -2.0, x_max = 2.0;
    double y_min = -2.0, y_max = 2.0;
    
    // Create coordinate vectors
    Eigen::VectorXd xx = Eigen::VectorXd::LinSpaced(nx, x_min, x_max);
    Eigen::VectorXd yy = Eigen::VectorXd::LinSpaced(ny, y_min, y_max);
    
    // ========================================================================
    // DEFINE YOUR TRIANGLE HERE
    // ========================================================================
    Triangle2D tri(
        -1.0, -0.5,   // Vertex 1 (bottom-left)
         1.0, -0.5,   // Vertex 2 (bottom-right)
         0.0,  1.0    // Vertex 3 (top)
    );
    
    // Print triangle info
    std::cout << "\nTriangle vertices:" << std::endl;
    std::cout << "  P1: (" << tri.x1 << ", " << tri.y1 << ")" << std::endl;
    std::cout << "  P2: (" << tri.x2 << ", " << tri.y2 << ")" << std::endl;
    std::cout << "  P3: (" << tri.x3 << ", " << tri.y3 << ")" << std::endl;
    writeTriangleGmsh(tri, "triangle_distance.msh");
    // ========================================================================
    // SETUP EIKONAL SOLVER PARAMETERS
    // ========================================================================
    
    // For distance computation, we need:
    // - Isotropic medium: a = b = 1.0, c = 0.0
    // - Slowness (fun) = 1.0 (speed = 1.0)
    // - Source point: centroid of the triangle (or any point inside)
    
    Eigen::MatrixXd a(nx, ny), b(nx, ny), c(nx, ny), fun(nx, ny);
    a.setConstant(1.0);
    b.setConstant(1.0);
    c.setConstant(0.0);
    fun.setConstant(1.0);  // Slowness = 1.0
    
    // Source point (centroid of triangle)
    double x0 = (tri.x1 + tri.x2 + tri.x3) / 3.0;
    double y0 = (tri.y1 + tri.y2 + tri.y3) / 3.0;
    
    std::cout << "\nEikonal solver parameters:" << std::endl;
    std::cout << "  Source point: (" << x0 << ", " << y0 << ")" << std::endl;
    std::cout << "  Medium: isotropic (a=b=1, c=0)" << std::endl;
    std::cout << "  Slowness: 1.0" << std::endl;
    
    // ========================================================================
    // COMPUTE DISTANCE FIELD USING EIKONAL SOLVER
    // ========================================================================
    
    Eigen::MatrixXd T_solver(nx, ny);
    bool use_cuda = false;  // Set to true if CUDA is available
    
    std::cout << "\nRunning Eikonal solver (FSM_WENO3_PS_2d)..." << std::endl;
    std::cout << "Grid size: " << nx << " x " << ny << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // CALL THE EIKONAL SOLVER
    EikonSolver2D::FSM_WENO3_PS_2d(xx, yy, a, b, c, fun, x0, y0, T_solver, use_cuda);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "Solver completed in " << duration.count() << " ms" << std::endl;
    
    // ========================================================================
    // COMPUTE REFERENCE DISTANCE (for comparison)
    // ========================================================================
    
    Eigen::MatrixXd T_reference(nx, ny);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            T_reference(i, j) = distanceToTriangle(xx(i), yy(j), tri);
        }
    }
    
    // ========================================================================
    // DEBUG OUTPUT
    // ========================================================================
    
    // Find statistics
    double minVal_solver = T_solver.minCoeff();
    double maxVal_solver = T_solver.maxCoeff();
    double minVal_ref = T_reference.minCoeff();
    double maxVal_ref = T_reference.maxCoeff();
    
    std::cout << "\nSolver distance field statistics:" << std::endl;
    std::cout << "  Min distance: " << minVal_solver << std::endl;
    std::cout << "  Max distance: " << maxVal_solver << std::endl;
    
    std::cout << "\nReference distance field statistics:" << std::endl;
    std::cout << "  Min distance: " << minVal_ref << std::endl;
    std::cout << "  Max distance: " << maxVal_ref << std::endl;
    
    // Sample distances
    std::cout << "\nSample distances (Solver vs Reference):" << std::endl;
    
    double cx = (tri.x1 + tri.x2 + tri.x3) / 3.0;
    double cy = (tri.y1 + tri.y2 + tri.y3) / 3.0;
    int idx_cx = static_cast<int>((cx - xx(0)) / (xx(1) - xx(0)) + 0.5);
    int idx_cy = static_cast<int>((cy - yy(0)) / (yy(1) - yy(0)) + 0.5);
    std::cout << "  Center (" << cx << ", " << cy << "): " 
              << T_solver(idx_cx, idx_cy) << " (solver), "
              << T_reference(idx_cx, idx_cy) << " (reference)" << std::endl;
    
    // Check a point outside
    double px = 0.0, py = -1.5;
    int idx_px = static_cast<int>((px - xx(0)) / (xx(1) - xx(0)) + 0.5);
    int idx_py = static_cast<int>((py - yy(0)) / (yy(1) - yy(0)) + 0.5);
    std::cout << "  Outside (" << px << ", " << py << "): " 
              << T_solver(idx_px, idx_py) << " (solver), "
              << T_reference(idx_px, idx_py) << " (reference)" << std::endl;
    
    // ========================================================================
    // SAVE RESULTS
    // ========================================================================
    
    saveAsVTK(T_solver, xx, yy, "triangle_distance_solver.vtk", "Distance_Solver");
    saveAsVTK(T_reference, xx, yy, "triangle_distance_reference.vtk", "Distance_Reference");
    
    // Compute difference
    Eigen::MatrixXd T_diff = (T_solver - T_reference).cwiseAbs();
    saveAsVTK(T_diff, xx, yy, "triangle_distance_diff.vtk", "Difference");
    
    std::cout << "\n==========================================" << std::endl;
    std::cout << "Example completed successfully!" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << "Files saved:" << std::endl;
    std::cout << "  - triangle_distance_solver.vtk (Eikonal solver result)" << std::endl;
    std::cout << "  - triangle_distance_reference.vtk (Direct computation)" << std::endl;
    std::cout << "  - triangle_distance_diff.vtk (Difference)" << std::endl;
    
    return 0;
}

