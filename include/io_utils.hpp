#ifndef IO_UTILS_H
#define IO_UTILS_H

#include <string>
#include <Eigen/Dense>
#include "polygon_utils.hpp"

// Read OBJ file
void readInputPolygon(std::string obj_filename,Polygon2D& poly);
Polygon2D readOBJFile(const std::string& filename);

// Write grid
void writeGridToVTK(const Eigen::VectorXd& xx, const Eigen::VectorXd& yy,
                    const std::string& filename);
void writeGridToMSH(const Eigen::VectorXd& xx, const Eigen::VectorXd& yy,
                    const std::string& filename);

// Write solution
void saveAsVTK(const Eigen::MatrixXd& T, const Eigen::VectorXd& xx, 
               const Eigen::VectorXd& yy, const std::string& filename,
               const std::string& fieldname);

// Write polygon geometry
void writePolygonGmsh(const Polygon2D& poly, const std::string& filename);
void writePolygonVTK(const Polygon2D& poly, const std::string& filename);

#endif // IO_UTILS_H