// PrismMesh2D.h
#ifndef PRISM_MESH_2D_H
#define PRISM_MESH_2D_H

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <set>
#include <memory>
#include"polygon_utils.hpp"


struct EdgeHash {
    size_t operator()(const Edge& e) const {
        return std::hash<int>()(e.v1) ^ (std::hash<int>()(e.v2) << 1);
    }
};

class PrismMesh2D {
public:
    PrismMesh2D();
    ~PrismMesh2D();

    // Vertex management using Point2D
    int addVertex(const Point2D& point);
    void setVertex(int idx, const Point2D& point);
    Point2D getVertex(int idx) const;
    int getNumVertices() const;
    void clear();

    // Normal computation and propagation
    Point2D computeVertexNormal(int vertexIdx, 
                                const std::vector<int>& neighborIndices) const;
    Point2D propagateVertexAlongNormal(int vertexIdx, 
                                       double distance,
                                       const std::vector<int>& neighborIndices);
    
    // Mesh generation
    void generatePrismMesh(const std::vector<std::vector<int>>& baseTriangles,
                          const std::vector<double>& extrudeDistances,
                          double topOffset = 0.0);
    void createVerts(Polygon2D& poly, int numlayer, std::vector<double>& dis);
    
  
    // File I/O
    void writeMSH(const std::string& filename) const;
    void writeMSHAdvanced(const std::string& filename, 
                         const std::map<std::string, std::vector<int>>& physicalGroups = {}) const;

    template<typename T, typename = void>
    struct has_print : std::false_type {};

    template<typename T>
    struct has_print<T, std::void_t<decltype(std::declval<T>().print())>> 
        : std::true_type {};

    template<typename T>
    void printVector(const std::vector<T>& vec, 
                const std::string& name = "Vector");
    

    // Visualization and statistics
    void visualize() const;
    void visualizeDetailed() const;
    
    struct MeshStatistics {
        int numVertices;
        int numPrisms;
        int numBoundaryEdges;
        std::vector<Point2D> vertexCoords;
        std::vector<std::vector<int>> prismConnectivity;
        std::vector<Edge> boundaryEdges;
    };
    MeshStatistics getMeshStatistics() const;

    // Getters
    const std::vector<Point2D>& getVertices() const { return vertices; }
    const std::vector<std::vector<int>>& getPrisms() const { return prisms; }
    const std::vector<Edge>& getBoundaryEdges() const { return boundaryEdges; }

private:
    std::vector<Point2D> vertices;
    std::map<int,std::vector<int>> v2v;
    std::map<int,std::vector<int>> v2e;
    std::vector<Edge> edges;
    std::map<int,std::vector<int>> e2e;
    std::map<int,std::vector<int>> e2f;
    

    std::vector<std::vector<int>> prisms;
    std::vector<Edge> boundaryEdges;
    mutable std::map<int, Point2D> vertexNormals; // Mutable for lazy computation

    // Helper functions
    std::vector<int> findVertexNeighbors(int vertexIdx, 
                                        const std::vector<std::vector<int>>& triangles) const;
    void addBoundaryEdges(const std::vector<int>& triangle);
    bool isEdgeOnBoundary(const Edge& edge, const std::vector<std::vector<int>>& triangles) const;
    bool isEdgeInPrism(const Edge& edge, const std::vector<int>& prism) const;
};

#endif // PRISM_MESH_2D_H