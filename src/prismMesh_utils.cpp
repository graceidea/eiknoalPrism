// PrismMesh2D.cpp
#include "prismMesh_utils.hpp"
#include <sstream>
#include <iomanip>
#include <type_traits>
PrismMesh2D::PrismMesh2D() {
    // Constructor
}

PrismMesh2D::~PrismMesh2D() {
    // Destructor
}

void PrismMesh2D::clear() {
    vertices.clear();
    prisms.clear();
    boundaryEdges.clear();
    vertexNormals.clear();
}
void PrismMesh2D::createVerts(Polygon2D& poly, int numlayer, std::vector<double>& dis)
{
    std::cout<<"\n======================================"<<std::endl;
    std::cout<<"    4. GENERATE PRISM MESH"<<std::endl;
    std::cout<<"======================================"<<std::endl;
    bool test=false;
    vertices.clear();
    if (numlayer<1 || dis.empty()) {
        std::cout<<" WARNING::numlayer="<<numlayer
            <<" disFunc="<<dis.size()<<std::endl;
        return;
    }
    if (static_cast<int>(dis.size())!=numlayer) {
        std::cout<<"dis.size="<<dis.size()<<" !=numlayer! numlayey="
           <<numlayer<<std::endl;
        return;
    }
    Vector2D norm;

    for(unsigned j=0;j<static_cast<unsigned>(numlayer);++j) {
      Point2D newp;
      for (unsigned i=0;i<poly.vertices.size();++i) {
        norm=poly.compute_node_normal(i);
        if (test) norm.print();

        v2v[static_cast<int>(vertices.size())].push_back(numlayer);
        Point2D pt = getVertex(i);
        vertices.push_back(pt);
     
        Point2D newPT;
        newPT.x=pt.x+norm.x*dis[j];
        newPT.y=pt.y+norm.y*dis[j];
        newPT.setId(vertices.size());
        vertices.push_back(newPT);
        poly.setVertex(newPT,i);//reset poly vert
        Edge edge(edges.size(),vertices.size()-2,vertices.size()-1);
        edges.push_back(edge);
     }     
    }
    this->printVector(vertices,"vertices");
    this->printVector(edges,"edges");
}
int PrismMesh2D::addVertex(const Point2D& point) {
    Point2D newPoint = point;
    newPoint.id = static_cast<int>(vertices.size());
    vertices.push_back(newPoint);
    return static_cast<int>(vertices.size()) - 1;
}

void PrismMesh2D::setVertex(int idx, const Point2D& point) {
    if (idx >= 0 && idx < static_cast<int>(vertices.size())) {
        vertices[idx] = point;
        vertices[idx].id = idx;
    }
}

Point2D PrismMesh2D::getVertex(int idx) const {
    if (idx >= 0 && idx < static_cast<int>(vertices.size())) {
        return vertices[idx];
    }
    return Point2D(0.0, 0.0, -1);
}

int PrismMesh2D::getNumVertices() const {
    return vertices.size();
}

Point2D PrismMesh2D::computeVertexNormal(int vertexIdx, 
                                         const std::vector<int>& neighborIndices) const {
    if (neighborIndices.size() < 2) {
        return Point2D(0.0, 0.0);
    }
    
    Point2D normalSum(0.0, 0.0);
    
    for (size_t i = 0; i < neighborIndices.size(); ++i) {
        size_t j = (i + 1) % neighborIndices.size();
        
        Point2D v1 = vertices[neighborIndices[i]];
        Point2D v2 = vertices[neighborIndices[j]];
        
        // Edge vector
        Point2D edge = v2 - v1;
        
        // Rotate 90 degrees to get normal
        Point2D normal(-edge.y, edge.x);
        
        // Ensure normal points outward
        Point2D center = vertices[vertexIdx];
        Point2D midPoint = (v1 + v2) * 0.5;
        Point2D toCenter = center - midPoint;
        
        if (normal.dot(toCenter) < 0) {
            normal = normal * (-1.0);
        }
        
        // Normalize and add to sum
        normal = normal.normalized();
        normalSum = normalSum + normal;
    }
    
    // Average and normalize
    return normalSum.normalized();
}

Point2D PrismMesh2D::propagateVertexAlongNormal(int vertexIdx, 
                                                double distance,
                                                const std::vector<int>& neighborIndices) {
    // Compute normal at the vertex
    Point2D normal = computeVertexNormal(vertexIdx, neighborIndices);
    
    if (normal.norm() < 1e-10) {
        throw std::runtime_error("Could not compute normal for vertex " + std::to_string(vertexIdx));
    }
    
    // Propagate the vertex
    Point2D currentPos = vertices[vertexIdx];
    Point2D newPos = currentPos + (normal * distance);
    
    // Update the vertex position
    vertices[vertexIdx] = newPos;
    vertices[vertexIdx].id = vertexIdx;
    
    return newPos;
}

std::vector<int> PrismMesh2D::findVertexNeighbors(int vertexIdx, 
                                                  const std::vector<std::vector<int>>& triangles) const {
    std::set<int> neighborSet;
    
    for (const auto& triangle : triangles) {
        if (std::find(triangle.begin(), triangle.end(), vertexIdx) != triangle.end()) {
            for (int v : triangle) {
                if (v != vertexIdx) {
                    neighborSet.insert(v);
                }
            }
        }
    }
    
    return std::vector<int>(neighborSet.begin(), neighborSet.end());
}

bool PrismMesh2D::isEdgeOnBoundary(const Edge& edge, 
                                   const std::vector<std::vector<int>>& triangles) const {
    int edgeCount = 0;
    
    for (const auto& triangle : triangles) {
        for (int i = 0; i < 3; ++i) {
            int v1 = triangle[i];
            int v2 = triangle[(i + 1) % 3];
            if ((v1 == edge.v1 && v2 == edge.v2) || (v1 == edge.v2 && v2 == edge.v1)) {
                edgeCount++;
                if (edgeCount > 1) {
                    return false;
                }
                break;
            }
        }
    }
    
    return edgeCount == 1;
}

bool PrismMesh2D::isEdgeInPrism(const Edge& edge, const std::vector<int>& prism) const {
    // Check if edge is in the bottom or top face of the prism
    for (int i = 0; i < 3; ++i) {
        int v1 = prism[i];
        int v2 = prism[(i + 1) % 3];
        if ((v1 == edge.v1 && v2 == edge.v2) || (v1 == edge.v2 && v2 == edge.v1)) {
            return true;
        }
    }
    return false;
}

void PrismMesh2D::addBoundaryEdges(const std::vector<int>& triangle) {
    std::vector<Edge> edges = {
        Edge(triangle[0], triangle[1]),
        Edge(triangle[1], triangle[2]),
        Edge(triangle[2], triangle[0])
    };
    
    // Convert prisms to triangles for edge checking
    std::vector<std::vector<int>> baseTriangles;
    for (const auto& prism : prisms) {
        std::vector<int> bottomFace = {prism[0], prism[1], prism[2]};
        baseTriangles.push_back(bottomFace);
    }
    
    for (const auto& edge : edges) {
        if (isEdgeOnBoundary(edge, baseTriangles)) {
            // Check if this edge already exists in boundaryEdges
            bool exists = false;
            for (const auto& be : boundaryEdges) {
                if (be == edge) {
                    exists = true;
                    break;
                }
            }
            if (!exists) {
                boundaryEdges.push_back(edge);
            }
        }
    }
}

void PrismMesh2D::generatePrismMesh(const std::vector<std::vector<int>>& baseTriangles,
                                   const std::vector<double>& extrudeDistances,
                                   double topOffset) {
    // Clear existing mesh
    prisms.clear();
    boundaryEdges.clear();
    vertexNormals.clear();
    
    int numBaseVertices = vertices.size();
    
    // Create top vertices by extruding base vertices
    for (int i = 0; i < numBaseVertices; ++i) {
        double dist = (i < static_cast<int>(extrudeDistances.size())) ? 
                       extrudeDistances[i] : extrudeDistances[0];
        
        // Find neighbors for this vertex
        std::vector<int> neighbors = findVertexNeighbors(i, baseTriangles);
        Point2D normal = computeVertexNormal(i, neighbors);
        
        if (normal.norm() < 1e-10) {
            // If normal cannot be computed, use vertical extrusion
            normal = Point2D(0.0, 1.0);
        }
        
        // Create top vertex
        Point2D topVertex = vertices[i] + (normal * (dist + topOffset));
        topVertex.id = vertices.size(); // Set id before adding
        addVertex(topVertex);
    }
    
    // Create prism elements
    for (const auto& triangle : baseTriangles) {
        if (triangle.size() != 3) {
            throw std::runtime_error("Base triangles must have exactly 3 vertices");
        }
        
        int b0 = triangle[0];
        int b1 = triangle[1];
        int b2 = triangle[2];
        
        // Top face vertices (newly created)
        int t0 = b0 + numBaseVertices;
        int t1 = b1 + numBaseVertices;
        int t2 = b2 + numBaseVertices;
        
        // Add prism (hexahedron with 6 vertices)
        // Order: bottom face CCW, top face CCW
        std::vector<int> prism = {b0, b1, b2, t0, t1, t2};
        prisms.push_back(prism);
        
        // Add boundary edges if this is an exterior face
        addBoundaryEdges(triangle);
    }
}

template<typename T>
void PrismMesh2D::printVector(const std::vector<T>& vec, 
    const std::string& name)
{
    std::cout << name << " (size=" << vec.size() << "):\n";
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << "  [" << i << "] ";
        if constexpr (has_print<T>::value) {
            vec[i].print();
        } else {
            std::cout << vec[i] << "\n";
        }
    }
}
void PrismMesh2D::writeMSH(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    
    // Write header
    file << "$MeshFormat\n";
    file << "2.2 0 8\n";  // Version 2.2, ASCII, 8-byte size
    file << "$EndMeshFormat\n";
    
    // Write vertices
    file << "$Nodes\n";
    file << vertices.size() << "\n";
    for (size_t i = 0; i < vertices.size(); ++i) {
        file << std::fixed << std::setprecision(10) 
             << (i + 1) << " " << vertices[i].x << " " << vertices[i].y << " 0.0\n";
    }
    file << "$EndNodes\n";
    
    // Write elements
    file << "$Elements\n";
    int totalElements = prisms.size() + boundaryEdges.size();
    file << totalElements << "\n";
    
    int elementId = 1;
    // Write prism elements
    for (const auto& prism : prisms) {
        file << elementId << " 6 2 1 1";  // Element type 6 = hexahedron
        for (int v : prism) {
            file << " " << (v + 1);
        }
        file << "\n";
        elementId++;
    }
    
    // Write boundary edges
    for (const auto& edge : boundaryEdges) {
        file << elementId << " 1 2 1 1";  // Element type 1 = line
        file << " " << (edge.v1 + 1) << " " << (edge.v2 + 1);
        file << "\n";
        elementId++;
    }
    
    file << "$EndElements\n";
    file.close();
}

void PrismMesh2D::writeMSHAdvanced(const std::string& filename, 
                                   const std::map<std::string, std::vector<int>>& physicalGroups) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    
    // Write header
    file << "$MeshFormat\n";
    file << "2.2 0 8\n";
    file << "$EndMeshFormat\n";
    
    // Write vertices
    file << "$Nodes\n";
    file << vertices.size() << "\n";
    for (size_t i = 0; i < vertices.size(); ++i) {
        file << std::fixed << std::setprecision(10) 
             << (i + 1) << " " << vertices[i].x << " " << vertices[i].y << " 0.0\n";
    }
    file << "$EndNodes\n";
    
    // Write elements
    file << "$Elements\n";
    int totalElements = prisms.size() + boundaryEdges.size();
    file << totalElements << "\n";
    
    int elementId = 1;
    // Write prism elements
    for (const auto& prism : prisms) {
        file << elementId << " 6 2 1 1";
        for (int v : prism) {
            file << " " << (v + 1);
        }
        file << "\n";
        elementId++;
    }
    
    // Write boundary edges
    for (const auto& edge : boundaryEdges) {
        file << elementId << " 1 2 1 1";
        file << " " << (edge.v1 + 1) << " " << (edge.v2 + 1);
        file << "\n";
        elementId++;
    }
    
    file << "$EndElements\n";
    
    // Write physical groups if provided
    if (!physicalGroups.empty()) {
        file << "$PhysicalNames\n";
        file << physicalGroups.size() << "\n";
        for (const auto& group : physicalGroups) {
            file << "1 " << group.second.size() << " \"" << group.first << "\"\n";
        }
        file << "$EndPhysicalNames\n";
    }
    
    file.close();
}

PrismMesh2D::MeshStatistics PrismMesh2D::getMeshStatistics() const {
    MeshStatistics stats;
    stats.numVertices = vertices.size();
    stats.numPrisms = prisms.size();
    stats.numBoundaryEdges = boundaryEdges.size();
    stats.vertexCoords = vertices;
    stats.prismConnectivity = prisms;
    stats.boundaryEdges = boundaryEdges;
    return stats;
}

void PrismMesh2D::visualize() const {
    std::cout << "\n=== Prism Mesh Visualization (ASCII) ===\n";
    std::cout << "Vertices:\n";
    for (size_t i = 0; i < vertices.size(); ++i) {
        std::cout << "  " << i << ": (" << vertices[i].x << ", " << vertices[i].y 
                  << ") [id=" << vertices[i].id << "]\n";
    }
    
    std::cout << "\nPrisms:\n";
    for (size_t i = 0; i < prisms.size(); ++i) {
        std::cout << "  " << i << ": [";
        for (size_t j = 0; j < prisms[i].size(); ++j) {
            std::cout << prisms[i][j];
            if (j < prisms[i].size() - 1) std::cout << ", ";
        }
        std::cout << "]\n";
    }
    
    std::cout << "\nBoundary Edges:\n";
    for (size_t i = 0; i < boundaryEdges.size(); ++i) {
        std::cout << "  " << i << ": (" << boundaryEdges[i].v1 << ", " << boundaryEdges[i].v2 << ")\n";
    }
}

void PrismMesh2D::visualizeDetailed() const {
    visualize();
    
    // Additional detailed statistics
    std::cout << "\n=== Detailed Mesh Information ===\n";
    std::cout << "Total vertices: " << vertices.size() << "\n";
    std::cout << "Total prisms: " << prisms.size() << "\n";
    std::cout << "Total boundary edges: " << boundaryEdges.size() << "\n";
    
    // Calculate average prism aspect ratio (simple measure)
    double avgAspect = 0.0;
    for (const auto& prism : prisms) {
        Point2D p0 = vertices[prism[0]];
        Point2D p1 = vertices[prism[1]];
        Point2D p2 = vertices[prism[2]];
        //Point2D p3 = vertices[prism[3]];
        
        double area = std::abs((p1 - p0).cross(p2 - p0)) / 2.0;
        if (area > 0) {
            avgAspect += area;
        }
    }
    if (!prisms.empty()) {
        avgAspect /= prisms.size();
        std::cout << "Average prism area: " << avgAspect << "\n";
    }
}