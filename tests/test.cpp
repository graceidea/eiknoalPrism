#include <iostream>
#include <string>
#include <chrono>
#include <fstream>
#include <limits>
#include <cmath>
#include <Eigen/Dense>

#include "eikon_solver_2d.hpp"
#include "polygon_utils.hpp"
#include "mesh_utils.hpp"
#include "io_utils.hpp"


int main(int argc, char** argv) {
    std::string obj_filename = "input.obj";
    if (argc > 1) {
        obj_filename = argv[1];
    }
    Polygon2D poly;
    readInputPolygon(obj_filename,poly);//,grid);
  
    //int nx, ny;
    //double xmin, xmax, ymin, ymax;
    //double grid_spacing = 0.05;  // Adjust as needed
    Grid grid;
    computeGridParameters(poly,grid);// nx, ny, xmin, xmax, ymin, ymax, grid_spacing);
    
    Eigen::MatrixXd T_exact(grid.nx, grid.ny);
    Eigen::VectorXd xx = Eigen::VectorXd::LinSpaced(grid.nx, grid.xmin, grid.xmax);
    Eigen::VectorXd yy = Eigen::VectorXd::LinSpaced(grid.ny, grid.ymin, grid.ymax);         
    computeDistance(poly,T_exact,xx,yy,grid);
   
    writeOutput(poly, grid, T_exact,xx, yy);   
    
    return 0;
}
