#ifndef POLYGON_UTILS_H
#define POLYGON_UTILS_H

#include <vector>
#include <cmath>
#include <iostream>
struct Point2D {
    double x, y;
    int id;
    Point2D(double x_ = 0.0, double y_ = 0.0, int id_=1) 
     : x(x_), y(y_), id(id_) {}
    void setId(int id_) {id=id_;}
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
    Point2D centroid() const;
    void boundingBox(double& xmin, double& xmax, double& ymin, double& ymax) const;
    double area() const;
    double distanceToBoundary(const Point2D& p) const;
    bool isOnBoundary(const Point2D& p, double tol) const;
    bool contains(const Point2D& p) const;
    bool isValid() const;
    size_t numEdges() const { return vertices.size(); }
    double perimeter() const;   

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