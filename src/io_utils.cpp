#include "io_utils.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <set>
#include <map>

void readInputPolygon(std::string obj_filename,Polygon2D& poly)
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

Polygon2D readOBJFile(const std::string& filename) 
{
    bool test=false;
    Polygon2D poly;
    std::ifstream file(filename);
    std::cout<<"  Polygon2D::readOBJFile"<<std::endl;
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file: " << filename << std::endl;
        return poly;
    }
    
    std::string line;
    std::vector<Point2D> all_points;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;
        
        if (prefix == "v") {
            double x, y, z = 0.0;
            iss >> x >> y >> z;
            all_points.push_back(Point2D(x, y));
        }
        else if (prefix == "f") {           
            Face f;
            f.vertices.clear();
            std::string token;
            while (iss >> token) {
                std::string v_str = token;
                size_t pos = v_str.find('/');
                if (pos != std::string::npos) {
                    v_str = v_str.substr(0, pos);
                }
                if (!v_str.empty()) {
                    int idx = std::stoi(v_str) - 1;
                    if (idx >= 0 && idx < static_cast<int>(all_points.size())) 
                    {
                      if (test) std::cout<<" idx="<<idx;
                      f.vertices.push_back(idx);
                    }
                    //else if (test) std::cout<<"  f.vertices.size="<<f.vertices.size();
                }
            }
            test=true;
            if (test)std::cout<<std::endl;
            if (!f.vertices.empty()) {
                if (test) {
                  std::cout<<" faceVert1="<<f.vertices.size()<<"(";
                  for (size_t ii = 0; ii < f.vertices.size(); ++ii)
                    std::cout<<" "<<f.vertices[ii];
                  std::cout<<")"<<std::endl;
                }
                poly.faces.push_back(f);

                // Generate edges for this face and add to poly.edges
                for (size_t i = 0; i < f.vertices.size(); ++i) {
                    int v1 = f.vertices[i];
                    int v2 = f.vertices[(i + 1) % f.vertices.size()]; // Wrap around to close the polygon
                    
                    // Create edge with unique ID (you can use poly.edges.size() as ID)
                    Edge newEdge(v1, v2, poly.edges.size());
                    
                    // Check if this edge already exists (avoid duplicates)
                    bool edgeExists = false;
                    for (const auto& existingEdge : poly.edges) {
                        if (existingEdge == newEdge) {
                            edgeExists = true;
                            break;
                        }
                    }
                    
                    if (!edgeExists) {
                        poly.edges.push_back(newEdge);
                        if (test) {
                            std::cout << "Added edge: ";
                            newEdge.print();
                        }
                    }
                }
            }
            test=false;
        }
    }
    file.close();
    
    if (!poly.faces.empty() && !all_points.empty()) {
        std::set<int> usedVertices;
        for (const auto& face : poly.faces) {
            for (int idx : face.vertices) {
                if (idx >= 0 && idx < static_cast<int>(all_points.size())) {
                    usedVertices.insert(idx);
                }
            }
        }
        // Add vertices to polygon
        for (int idx : usedVertices) {
            Point2D p = all_points[idx];
            p.id = static_cast<int>(poly.vertices.size());
            poly.vertices.push_back(p);
        }
        // Store faces with remapped vertex indices
        std::map<int, int> indexMap;
        for (size_t i = 0; i < poly.vertices.size(); ++i) {
            // Need to map original indices to new indices
            for (const auto& face : poly.faces) {
                for (int idx : face.vertices) {
                    if (idx == static_cast<int>(i)) {
                        indexMap[i] = static_cast<int>(poly.vertices.size()) - 1;
                    }
                }
            }
        }
         
        // Extract edges from faces
        //poly.extractEdgesFromFaces();
        
        std::cout << "  Loaded " << all_points.size() << " vertices, "
                  <<poly.vertices.size()<<" polyVerts, "
                  << poly.faces.size() << " faces, face[0] has "
                  <<poly.faces[0].vertices.size()<<" Verts "
                  << poly.edges.size() << " edges\n";    
    }
    return poly;
}

void writeGridToVTK(const Eigen::VectorXd& xx, const Eigen::VectorXd& yy,
                    const std::string& filename) 
{  
    bool test=false;
    int nx = xx.size();
    int ny = yy.size();
    
    std::ofstream vtk_file(filename);
    if (!vtk_file.is_open()) {
        std::cerr << "Error: Could not open " << filename << std::endl;
        return;
    }
    
    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "2D Grid\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET STRUCTURED_GRID\n";
    vtk_file << "DIMENSIONS " << nx << " " << ny << " 1\n";
    vtk_file << "POINTS " << nx * ny << " float\n";
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            vtk_file << xx(i) << " " << yy(j) << " 0.0\n";
        }
    }
    
    vtk_file.close();
    if (test)std::cout << "Grid saved: " << filename << std::endl;
}
void writeGridToMSH(const Eigen::VectorXd& xx, const Eigen::VectorXd& yy,
                        const std::string& filename) 
{  
    bool test=false;
    int nx = xx.size();
    int ny = yy.size();
    
    std::ofstream msh_file(filename);
    if (!msh_file.is_open()) {
        std::cerr << "Error: Could not open " << filename << std::endl;
        return;
    }
    
    // Set precision for output
    msh_file << std::setprecision(15) << std::scientific;
    
    // ============================================
    // HEADER
    // ============================================
    msh_file << "$MeshFormat\n";
    msh_file << "2.2 0 8\n";
    msh_file << "$EndMeshFormat\n";
    
    // ============================================
    // NODES
    // ============================================
    msh_file << "$Nodes\n";
    int num_nodes = nx * ny;
    msh_file << num_nodes << "\n";
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int node_id = j * nx + i + 1;
            msh_file << node_id << " " 
                     << xx(i) << " " 
                     << yy(j) << " 0.0\n";
        }
    }
    msh_file << "$EndNodes\n";
    
    // ============================================
    // ELEMENTS (quadrilaterals)
    // ============================================
    msh_file << "$Elements\n";
    int num_quads = (nx - 1) * (ny - 1);
    msh_file << num_quads << "\n";
    
    int elem_id = 1;
    for (int j = 0; j < ny - 1; ++j) {
        for (int i = 0; i < nx - 1; ++i) {
            // Node indices (1-based)
            int n1 = j * nx + i + 1;        // bottom-left
            int n2 = j * nx + (i + 1) + 1;  // bottom-right
            int n3 = (j + 1) * nx + (i + 1) + 1; // top-right
            int n4 = (j + 1) * nx + i + 1;  // top-left
            
            // Format: elem_id element_type num_tags tag1 tag2 n1 n2 n3 n4
            // element_type 3 = quadrilateral
            // num_tags 2 (physical tag, elementary tag)
            msh_file << elem_id << " 3 2 1 1 " 
                     << n1 << " " << n2 << " " << n3 << " " << n4 << "\n";
            elem_id++;
        }
    }
    msh_file << "$EndElements\n";
    
    msh_file.close();
    if (test)std::cout << "Grid saved: " << filename << std::endl;
}
void writeGridToMSH_22(const Eigen::VectorXd& xx, const Eigen::VectorXd& yy,
                        const std::string& filename) 
{  
    bool test = false;
    int nx = xx.size();
    int ny = yy.size();
    
    std::ofstream msh_file(filename);
    if (!msh_file.is_open()) {
        std::cerr << "Error: Could not open " << filename << std::endl;
        return;
    }
    
    // Gmsh MSH file format version 2.2 (simpler, widely supported)
    msh_file << "$MeshFormat\n";
    msh_file << "2.2 0 8\n";
    msh_file << "$EndMeshFormat\n";
    
    // ============================================
    // SECTION 1: NODES
    // ============================================
    msh_file << "$Nodes\n";
    int num_nodes = nx * ny;
    msh_file << num_nodes << "\n";
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int node_id = j * nx + i + 1;
            msh_file << node_id << " " << xx(i) << " " << yy(j) << " 0.0\n";
        }
    }
    msh_file << "$EndNodes\n";
    
    // ============================================
    // SECTION 2: ELEMENTS (quadrilaterals)
    // ============================================
    msh_file << "$Elements\n";
    int num_quads = (nx - 1) * (ny - 1);
    msh_file << num_quads << "\n";
    
    int elem_id = 1;
    for (int j = 0; j < ny - 1; ++j) {
        for (int i = 0; i < nx - 1; ++i) {
            int n1 = j * nx + i + 1;        // bottom-left
            int n2 = j * nx + (i + 1) + 1;  // bottom-right
            int n3 = (j + 1) * nx + (i + 1) + 1;  // top-right
            int n4 = (j + 1) * nx + i + 1;  // top-left
            
            // Format: elem_id type num_tags tag1 tag2 ... node1 node2 ...
            // Type 3 = quadrilateral, 2 tags (physical and elementary)
            msh_file << elem_id << " 3 2 1 " << elem_id << " " << n1 << " " << n2 << " " << n3 << " " << n4 << "\n";
            elem_id++;
        }
    }
    msh_file << "$EndElements\n";
    
    msh_file.close();
    if (test) std::cout << "Grid saved: " << filename << std::endl;
}
void saveAsVTK(const Eigen::MatrixXd& T, const Eigen::VectorXd& xx, 
               const Eigen::VectorXd& yy, const std::string& filename,
               const std::string& fieldname) 
{  
    bool test=false;
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
    if (test)std::cout << "Solution saved: " << filename << std::endl;
}

void writePolygonGmsh(const Polygon2D& poly, const std::string& filename) 
{
    bool test=false;
    std::ofstream msh_file(filename);
    if (!msh_file.is_open()) {
        std::cerr << "Error: Could not open " << filename << std::endl;
        return;
    }
    
    int n = poly.vertices.size();
    
    msh_file << "$MeshFormat\n";
    msh_file << "2.2 0 8\n";
    msh_file << "$EndMeshFormat\n";
    
    msh_file << "$Nodes\n";
    msh_file << n << "\n";
    for (int i = 0; i < n; i++) {
        msh_file << (i + 1) << " " << poly.vertices[i].x << " " << poly.vertices[i].y << " 0.0\n";
    }
    msh_file << "$EndNodes\n";
    
    msh_file << "$Elements\n";
    msh_file << n << "\n";
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        msh_file << (i + 1) << " 1 2 0 0 " << (i + 1) << " " << (j + 1) << "\n";
    }
    msh_file << "$EndElements\n";
    
    msh_file.close();
    if (test)std::cout << "Polygon saved (Gmsh): " << filename << std::endl;
}

void writePolygonVTK(const Polygon2D& poly, const std::string& filename) 
{
    bool test=false;
    std::ofstream vtk_file(filename);
    if (!vtk_file.is_open()) {
        std::cerr << "Error: Could not open " << filename << std::endl;
        return;
    }
    
    int n = poly.vertices.size();
    
    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "Polygon\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET POLYDATA\n";
    
    vtk_file << "POINTS " << n << " float\n";
    for (const auto& p : poly.vertices) {
        vtk_file << p.x << " " << p.y << " 0.0\n";
    }
    
    vtk_file << "LINES " << n << " " << (n * 3) << "\n";
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        vtk_file << "2 " << i << " " << j << "\n";
    }
    
    vtk_file.close();
    if (test)std::cout << "Polygon saved (VTK): " << filename << std::endl;
}