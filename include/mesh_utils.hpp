#ifndef MESH_UTILS_H
#define MESH_UTILS_H
#include <string>
#include <Eigen/Dense>
#include "polygon_utils.hpp"

struct Grid {
 int nx, ny;
 double xmin, xmax, ymin, ymax;
 double grid_spacing = 0.05; 
};

void readInputPolygon(std::string obj_filename,Polygon2D& poly);//,Grid& grid);
void computeGridParameters(const Polygon2D& poly, Grid& grid);
                           //int& nx, int& ny,
                           //double& xmin, double& xmax,
                           //double& ymin, double& ymax,
                           //double padding_factor = 0.2);
void computeDistance(const Polygon2D& poly,Eigen::MatrixXd& T_exact,
    Eigen::VectorXd& xx,Eigen::VectorXd& yy,Grid& grid);
void computeGradientMagnitude(Grid grid,Eigen::MatrixXd& T_exact,
     Eigen::VectorXd& xx, Eigen::VectorXd& yy,Eigen::MatrixXd& grad_mag);                           
//void ApplySolver(Polygon2D poly, Grid grid);
void writeOutput(Polygon2D poly, Grid grid, Eigen::MatrixXd& T_exact,
    Eigen::VectorXd& xx, Eigen::VectorXd& yy);
#endif // MESH_UTILS_H