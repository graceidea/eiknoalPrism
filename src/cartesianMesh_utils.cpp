#include "cartesianMesh_utils.hpp"
#include "io_utils.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>

// Check if OpenMP is available
#ifdef _OPENMP
    #include <omp.h>
    #define OMP_PARALLEL_FOR _Pragma("omp parallel for collapse(2)")
#else
    #define OMP_PARALLEL_FOR
#endif
void readInputPolygon(std::string obj_filename,Polygon2D& poly)//,Grid& grid)
{
    std::cout << "==========================================" << std::endl;
    std::cout << "  1. READ INPUT CLOSED CURVE FROM OBJ FILE" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    obj_filename = "input.obj";  
    
    std::cout << "  Reading polygon from: " << obj_filename << std::endl;
    // Check if file exists
    std::ifstream file_check(obj_filename);
    if (!file_check.good()) {
        std::cerr << "Error: Cannot open file " << obj_filename << std::endl;
        return;
    }
    file_check.close();
    poly = readOBJFile(obj_filename);
    
    if (poly.vertices.size() < 3) {
        std::cerr << "Error: Invalid polygon (need at least 3 vertices)" << std::endl;
        return;
    }
    //std::cout << "\n=== Parameters ===" << std::endl;
    std::cout << "  Polygon vertices: " << poly.vertices.size() << std::endl;
    std::cout << "  Polygon area: " << poly.area() << std::endl;    
}       

void computeGridParameters(const Polygon2D& poly, Grid& grid) 
{
    std::cout<<"\n===================================================="<<std::endl;
    std::cout<<"  2, COMPUTE ALL NECESSARY GRID PARAMETERS"<<std::endl;
    std::cout<<"===================================================="<<std::endl;
    
    poly.boundingBox(grid.xmin, grid.xmax, grid.ymin, grid.ymax);
    
    double dx = grid.xmax - grid.xmin;
    double dy = grid.ymax - grid.ymin;
    double pad_x = dx * grid.grid_spacing;
    double pad_y = dy * grid.grid_spacing;
    
    grid.xmin -= pad_x;
    grid.xmax += pad_x;
    grid.ymin -= pad_y;
    grid.ymax += pad_y;
    
    // Ensure square grid (use max extent)
    double extent = std::max(grid.xmax - grid.xmin, grid.ymax - grid.ymin);
    double center_x = (grid.xmin + grid.xmax) / 2.0;
    double center_y = (grid.ymin + grid.ymax) / 2.0;
    
    grid.xmin = center_x - extent / 2.0;
    grid.xmax = center_x + extent / 2.0;
    grid.ymin = center_y - extent / 2.0;
    grid.ymax = center_y + extent / 2.0;
    
    // Use grid resolution based on polygon complexity
    grid.nx = 201;
    grid.ny = 201;
    // Validate grid parameters
    if (grid.nx <= 0 || grid.ny <= 0) {
        throw std::runtime_error("Invalid grid dimensions");
    }
    std::cout << "  Grid: " << grid.nx << " x " << grid.ny << std::endl;
    std::cout << "  Domain: [" << grid.xmin << ", " << grid.xmax << "] x [" 
        << grid.ymin << ", " << grid.ymax << "]" << std::endl;
    std::cout << "  Grid spacing: " << grid.grid_spacing << std::endl;
}


void computeDistance(const Polygon2D& poly,Eigen::MatrixXd& T_exact,
    Eigen::VectorXd& xx,Eigen::VectorXd& yy,Grid& grid)
{
    std::cout<<"\n========================================================"<<std::endl;
    std::cout<<"  3. COMPUTE EXACT DISTANCE FIELD TO POLYGON BOUNDARY"<<std::endl;
    std::cout<<"     This IS the exact solution for |∇T| = 1 with T=0 on boundary"<<std::endl;
    std::cout<<"========================================================"<<std::endl;
    //std::cout << "\n=== Computing Exact Distance Field ===" << std::endl;
    std::cout << "  Computing distance to polygon boundary for all grid points..." << std::endl;
        
    auto start = std::chrono::high_resolution_clock::now(); 
        
    #ifdef _OPENMP
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < grid.nx; i++) {
          for (int j = 0; j < grid.ny; j++) {
            Point2D p(xx(i), yy(j));
            double dist = poly.distanceToBoundary(p);
            T_exact(i,j) = poly.contains(p) ? dist : -dist;
          }
        }
    #else
        // Sequential version without OpenMP
        for (int i = 0; i < grid.nx; i++) {
          for (int j = 0; j < grid.ny; j++) {
            Point2D p(xx(i), yy(j));
            double dist = poly.distanceToBoundary(p);
            T_exact(i,j) = poly.contains(p) ? dist : -dist;
          }
        }
    #endif

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "  Distance computation completed in " << duration.count() << " ms" << std::endl;
        
    // Compute some statistics
    double min_dist = T_exact.minCoeff();
    double max_dist = T_exact.maxCoeff();
    double mean_dist = T_exact.mean();
    std::cout << " Distance range: [" << min_dist << ", " << max_dist << "]" << std::endl;
    std::cout << " Mean distance: " << mean_dist << std::endl;
        
    // Verify boundary condition (points near boundary should have T ≈ 0)
    double boundary_error = 0.0;
    int boundary_count = 0;
    double tol = grid.grid_spacing * 0.5;
        
    for (int i = 0; i < grid.nx; i++) {
        for (int j = 0; j < grid.ny; j++) {
            Point2D p(xx(i), yy(j));
            if (poly.isOnBoundary(p, tol)) {
                boundary_error += std::abs(T_exact(i,j));
                boundary_count++;
            }
        }
    }
        
    if (boundary_count > 0) {
            boundary_error /= boundary_count;
            std::cout << " Average boundary T error: " << boundary_error << std::endl;
            if (boundary_error > 1e-6) {
                std::cout << " Warning: Boundary condition not accurately satisfied" << std::endl;
            }
    }
}
void computeGradientMagnitude(Grid grid,Eigen::MatrixXd& T_exact,
      Eigen::VectorXd& xx, Eigen::VectorXd& yy,Eigen::MatrixXd& grad_mag)
{
    // Compute and write gradient magnitude (should be ≈ 1)
    for (int i = 1; i < grid.nx-1; i++) {
        for (int j = 1; j < grid.ny-1; j++) {
            double dx = (T_exact(i+1,j) - T_exact(i-1,j)) / (xx(i+1) - xx(i-1));
            double dy = (T_exact(i,j+1) - T_exact(i,j-1)) / (yy(j+1) - yy(j-1));
            grad_mag(i,j) = std::sqrt(dx*dx + dy*dy);
        }
    }
    for (int i = 0; i < grid.nx; i++) {
        grad_mag(i,0) = grad_mag(i,1);
        grad_mag(i,grid.ny-1) = grad_mag(i,grid.ny-2);
    }
    for (int j = 0; j < grid.ny; j++) {
        grad_mag(0,j) = grad_mag(1,j);
        grad_mag(grid.nx-1,j) = grad_mag(grid.nx-2,j);
    }
}
void writeOutput(Polygon2D poly, Grid grid, Eigen::MatrixXd& T_exact,
    Eigen::VectorXd& xx, Eigen::VectorXd& yy)
{     
    std::cout<<"\n===================================="<<std::endl;                                                   
    std::cout << "  Writing Output Files " << std::endl;
    std::cout<<"===================================="<<std::endl; 

    // Write exact distance field
    saveAsVTK(T_exact, xx, yy, "distance_exact.vtk", "Distance");
    std::cout << " Wrote: distance_exact.vtk" << std::endl;
             
    Eigen::MatrixXd grad_mag(grid.nx, grid.ny); 
    computeGradientMagnitude(grid,T_exact,xx,yy,grad_mag);  
    saveAsVTK(grad_mag, xx, yy, "gradient_magnitude.vtk", "GradientMagnitude");
    std::cout << " Wrote: gradient_magnitude.vtk" << std::endl;

    writePolygonVTK(poly, "polygon.vtk");
    std::cout << " Wrote: polygon.vtk" << std::endl;       
    writePolygonGmsh(poly, "polygon.msh");
    std::cout << " Wrote: polygon.msh" << std::endl;
        
    writeGridToVTK(xx, yy, "grid.vtk");
    std::cout << " Wrote: grid.vtk" << std::endl;
    writeGridToMSH(xx, yy, "grid.msh");
    std::cout << " Wrote: grid.msh" << std::endl; 

    std::cout << "\n=== Results Summary ===" << std::endl;
    std::cout << " Exact solution is the distance to polygon boundary" << std::endl;
    std::cout << " |∇T| ≈ " << grad_mag.mean() << " (should be ≈ 1.0)" << std::endl;
        
    std::cout << "\n=== Output Files ===" << std::endl;
    std::cout << " 1. distance_exact.vtk      - Distance to polygon boundary" << std::endl;
    std::cout << " 2. gradient_magnitude.vtk  - |∇T| (should be ≈ 1)" << std::endl;
    std::cout << " 3. polygon.vtk             - Polygon boundary" << std::endl;
    std::cout << " 4. polygon.msh             - Polygon for Gmsh" << std::endl;
    std::cout << " 5. grid.vtk                - Grid points" << std::endl;
    std::cout << " 6. grid.msh                - Grid points" << std::endl;
        
    std::cout << "\n=== Visualization ===" << std::endl;
    std::cout << " ParaView: paraview distance_exact.vtk polygon.vtk" << std::endl;
    std::cout << " Gmsh:     gmsh distance_exact.vtk polygon.msh" << std::endl;
    std::cout << "\n==========================================" << std::endl;        
} 
