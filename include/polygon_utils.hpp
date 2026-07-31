#ifndef POLYGON_UTILS_H
#define POLYGON_UTILS_H

#include <vector>
#include <cmath>
#include <iostream>
struct Point2D {
    double x, y;
    int id;
    Point2D(double x_ = 0.0, double y_ = 0.0, int id_=-1) 
     : x(x_), y(y_), id(id_) {}
    void setId(int id_) {id=id_;}
    int getId() const {return id;}
    void print() const{ std::cout<<"  ptID="<<id<<"("<<x<<","<<y<<")"<<std::endl;}
    double distanceTo(const Point2D& other) const {
        double dx = x - other.x;
        double dy = y - other.y;
        return std::sqrt(dx*dx + dy*dy);
    }  
    double distanceSquaredTo(const Point2D& other) const {
        double dx = x - other.x;
        double dy = y - other.y;
        return dx*dx + dy*dy;
    }
    Point2D operator-(const Point2D& other) const {
        return Point2D(x - other.x, y - other.y);
    }   
    Point2D operator+(const Point2D& other) const {
        return Point2D(x + other.x, y + other.y);
    }
    Point2D operator*(double scalar) const {
        return Point2D(x * scalar, y * scalar);
    }    
    double cross(const Point2D& other) const {
        return x * other.y - y * other.x;
    }    
    double dot(const Point2D& other) const {
        return x * other.x + y * other.y;
    }
    double norm() const {
        return std::sqrt(x*x + y*y);
    }
    
    double normSquared() const {
        return x*x + y*y;
    }
    
    Point2D normalized() const {
        double n = norm();
        if (n > 1e-12) {
            return Point2D(x/n, y/n);
        }
        return Point2D(0, 0);
    }
   
};
// Define operator<< OUTSIDE the struct (as a free function)
inline std::ostream& operator<<(std::ostream& os, const Point2D& p) {
    os << "Point2D(id=" << p.id << ", x=" << p.x << ", y=" << p.y << ")";
    return os;
}
//=====================Edge
struct Edge {
    int v1, v2;     // Declared first
    int id;         // Declared second
    
    // Initializer list should match declaration order: v1, v2, id
    Edge(int v1_ = -1, int v2_ = -1, int i = -1) : v1(v1_), v2(v2_), id(i) {}
    
    bool operator==(const Edge& other) const {
        return (v1 == other.v1 && v2 == other.v2) || 
               (v1 == other.v2 && v2 == other.v1);
    }

    void print() const {
        std::cout << " EDGE=" << id << " edgeV1=" << v1 << " edgeV2=" << v2 << std::endl;
    }
};
struct Face {
    int v1, v2, v3;     // Declared first
    int e1, e2, e3;
    int id;         // Declared second
    std::vector<int> vertices;
    // Initializer list should match declaration order: v1, v2, id
    Face(int v1_ = -1, int v2_ = -1, int v3_=-1,int i = -1) 
      : v1(v1_), v2(v2_), v3(v3_), id(i) {
        vertices.push_back(v1);
        vertices.push_back(v2);
        vertices.push_back(v3);
      }
    
    void setFace(int id_, std::vector<int>& fv) {
        id=id_;
        for (unsigned i=0;i<fv.size();++i) {
         if (fv[i]>=0)     vertices.push_back(fv[i]);
         else std::cout<<"ERROR: faceid "<<fv[i]<<std::endl;
        }
    }
    void setFace(int v){ if (v>=0) vertices.push_back(v);}
    bool operator==(const Face& other) const {
        return (v1 == other.v1 && v2 == other.v2) || 
               (v1 == other.v2 && v2 == other.v1);
    }

    void print() const {
        std::cout << " FACE ID=" << id 
            <<" Verts=("<<v1<<","<<v2<<","<<v3<<")"
            <<" edges=("<<e1<<","<<e2<<","<<e3<<")"<< std::endl;
    }
};
// Free function operator<< for Edge
inline std::ostream& operator<<(std::ostream& os, const Edge& e) {
    os << "Edge(id=" << e.id << ", v1=" << e.v1 << ", v2=" << e.v2 << ")";
    return os;
}
//==================================end of Edge
// Structure representing a 2D Vector / Point
struct Vector2D {
    double x = 0.0;
    double y = 0.0;

    // Helper method to compute the magnitude (length) of the vector
    double length() const {
        return std::sqrt(x * x + y * y);
    }
    void print() const {
        std::cout << "(" << x << ", " << y << ")\n";
    }
    // Helper method to return a normalized (unit) vector
    Vector2D normalized() const {
        double len = length();
        if (len > 0.0) {
            return { x / len, y / len };
        }
        return { 0.0, 0.0 }; // Return zero vector if length is 0
    }
};
struct Polygon2D {
    std::vector<Point2D> vertices;
    std::vector<Edge> edges;
    std::vector<Face> faces;    
    Point2D centroid() const;
    Point2D getVertex(int idx) const;
    void setVertex(Point2D p, int i) 
    {
        if (i>=0&& i<static_cast<int>(vertices.size())) 
          vertices[i]=p;
    }
    void boundingBox(double& xmin, double& xmax, double& ymin, double& ymax) const;
    double area() const;
    double distanceToBoundary(const Point2D& p) const;
    bool isOnBoundary(const Point2D& p, double tol) const;
    bool contains(const Point2D& p) const;
    bool isValid() const;
    size_t numEdges() const { return vertices.size(); }
    double perimeter() const;   
    //void extractEdgesFromFaces();
    Vector2D compute_node_normal(size_t index) const;
private:
    // Helper for point-to-segment distance
    double pointToSegmentDistance(const Point2D& p, const Point2D& a, const Point2D& b) const;
    
    // Helper for point-on-segment test
    bool pointOnSegment(const Point2D& p, const Point2D& a, const Point2D& b, double tol) const;
    
    // Helper for segment intersection
    bool segmentsIntersect(const Point2D& a1, const Point2D& b1, 
                          const Point2D& a2, const Point2D& b2) const;     
};

// Distance functions
double pointToSegmentDistance(double px, double py, 
                              double x1, double y1, 
                              double x2, double y2);

bool pointInPolygon(double px, double py, const Polygon2D& poly);

double signedDistanceToPolygon(double px, double py, const Polygon2D& poly);

#endif // POLYGON_UTILS_H