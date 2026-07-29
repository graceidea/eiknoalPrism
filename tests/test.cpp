#include <iostream>
#include <string>
#include <chrono>
#include <fstream>
#include <limits>
#include <cmath>
#include <Eigen/Dense>

#include "eikon_solver_2d.hpp"
#include "polygon_utils.hpp"
#include "cartesianMesh_utils.hpp"
#include "io_utils.hpp"


int main(int argc, char** argv) {
    //1, readin input polygon as initial source
    std::string obj_filename = "input.obj";
    if (argc > 1) {
        obj_filename = argv[1];
    }
    Polygon2D poly;
    readInputPolygon(obj_filename,poly);//,grid);
  
    //2, compute background cartisian grid
    Grid grid;
    computeGridParameters(poly,grid);// nx, ny, xmin, xmax, ymin, ymax, grid_spacing);
    
    //3, compute signed distance 
    Eigen::MatrixXd T_exact(grid.nx, grid.ny);
    Eigen::VectorXd xx = Eigen::VectorXd::LinSpaced(grid.nx, grid.xmin, grid.xmax);
    Eigen::VectorXd yy = Eigen::VectorXd::LinSpaced(grid.ny, grid.ymin, grid.ymax);         
    computeDistance(poly,T_exact,xx,yy,grid);
   
    //3, write outputs
    writeOutput(poly, grid, T_exact,xx, yy);   
    Vector2D norm;
    for (unsigned i=0;i<poly.vertices.size();++i) {
     norm=poly.compute_node_normal(i);
     norm.print();
    }
    return 0;
}
