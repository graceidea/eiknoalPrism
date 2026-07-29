#include "polygon_utils.hpp"
#include <algorithm>
#include <limits>

Point2D Polygon2D::centroid() const {
    if (vertices.empty()) return Point2D(0, 0);
    double cx = 0.0, cy = 0.0;
    for (const auto& p : vertices) {
        cx += p.x;
        cy += p.y;
    }
    return Point2D(cx / vertices.size(), cy / vertices.size());
}

void Polygon2D::boundingBox(double& xmin, double& xmax, 
                            double& ymin, double& ymax) const {
    if (vertices.empty()) {
        xmin = xmax = ymin = ymax = 0.0;
        return;
    }
    xmin = xmax = vertices[0].x;
    ymin = ymax = vertices[0].y;
    for (const auto& p : vertices) {
        xmin = std::min(xmin, p.x);
        xmax = std::max(xmax, p.x);
        ymin = std::min(ymin, p.y);
        ymax = std::max(ymax, p.y);
    }
}

double Polygon2D::area() const {
    if (vertices.size() < 3) return 0.0;
    double a = 0.0;
    int n = vertices.size();
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        a += vertices[i].x * vertices[j].y;
        a -= vertices[j].x * vertices[i].y;
    }
    return 0.5 * a;
}

double Polygon2D::distanceToBoundary(const Point2D& p) const {
    if (vertices.size() < 3) {
        return std::numeric_limits<double>::infinity();
    }
    
    double minDist = std::numeric_limits<double>::infinity();
    
    // Compute distance to each edge
    for (size_t i = 0; i < vertices.size(); i++) {
        size_t j = (i + 1) % vertices.size();
        double dist = pointToSegmentDistance(p, vertices[i], vertices[j]);
        if (dist < minDist) {
            minDist = dist;
        }
    }
    return minDist;
}
double Polygon2D::pointToSegmentDistance(const Point2D& p, const Point2D& a, const Point2D& b) const {
    // Vector from a to b
    Point2D ab = b - a;
    Point2D ap = p - a;
    
    // Projection of ap onto ab, clamped to [0,1]
    double abLenSq = ab.normSquared();
    
    if (abLenSq < 1e-12) {
        // Segment is actually a point
        return p.distanceTo(a);
    }
    double t = ap.dot(ab) / abLenSq;
    t = std::max(0.0, std::min(1.0, t));
    
    // Closest point on segment
    Point2D closest = a + ab * t;
    
    return p.distanceTo(closest);
}
bool Polygon2D::isOnBoundary(const Point2D& p, double tol) const {
    if (vertices.size() < 3) {
        return false;
    }
    
    // Check if point lies on any edge
    for (size_t i = 0; i < vertices.size(); i++) {
        size_t j = (i + 1) % vertices.size();
        if (pointOnSegment(p, vertices[i], vertices[j], tol)) {
            return true;
        }
    }
    return false;
}
bool Polygon2D::pointOnSegment(const Point2D& p, const Point2D& a, 
                               const Point2D& b, double tol) const {
    // First check if point is within bounding box (quick rejection)
    double minX = std::min(a.x, b.x) - tol;
    double maxX = std::max(a.x, b.x) + tol;
    double minY = std::min(a.y, b.y) - tol;
    double maxY = std::max(a.y, b.y) + tol;
    
    if (p.x < minX || p.x > maxX || p.y < minY || p.y > maxY) {
        return false;
    }
    Point2D ab = b - a;
    Point2D ap = p - a;
    double cross = ab.cross(ap);
    
    if (std::abs(cross) > tol) {
        return false;
    }
    
    // Check if point is between a and b using dot product
    double dot = ab.dot(ap);
    if (dot < -tol || dot > ab.normSquared() + tol) {
        return false;
    }
    return true;
}
bool Polygon2D::contains(const Point2D& p) const {
    if (vertices.size() < 3) {
        return false;
    }
    
    bool inside = false;
    
    // Ray casting algorithm
    for (size_t i = 0, j = vertices.size() - 1; i < vertices.size(); j = i++) {
        const Point2D& vi = vertices[i];
        const Point2D& vj = vertices[j];
        bool intersect = ((vi.y > p.y) != (vj.y > p.y)) &&
                         (p.x < (vj.x - vi.x) * (p.y - vi.y) / (vj.y - vi.y) + vi.x);
        
        if (intersect) {
            inside = !inside;
        }
    }
    
    return inside;
}
bool Polygon2D::isValid() const {
    // Need at least 3 vertices
    if (vertices.size() < 3) {
        return false;
    }
    
    // Check for duplicate consecutive vertices
    for (size_t i = 0; i < vertices.size(); i++) {
        size_t j = (i + 1) % vertices.size();
        if (vertices[i].distanceTo(vertices[j]) < 1e-10) {
            return false;  // Duplicate vertex
        }
    }
   for (size_t i = 0; i < vertices.size(); i++) {
        size_t j = (i + 1) % vertices.size();
        for (size_t k = i + 1; k < vertices.size(); k++) {
            size_t l = (k + 1) % vertices.size();
            
            // Skip adjacent edges
            if (j == k || i == l) {
                continue;
            }
            if (segmentsIntersect(vertices[i], vertices[j], vertices[k], vertices[l])) {
                // Check if they just touch at endpoints (allowed)
                bool shareEndpoint = (vertices[i].distanceTo(vertices[k]) < 1e-10) ||
                                    (vertices[i].distanceTo(vertices[l]) < 1e-10) ||
                                    (vertices[j].distanceTo(vertices[k]) < 1e-10) ||
                                    (vertices[j].distanceTo(vertices[l]) < 1e-10);
                
                if (!shareEndpoint) {
                    return false;  // Self-intersection
                }
            }
        }
    }
    if (std::abs(area()) < 1e-12) {
        return false;
    }
    
    return true;
}
bool Polygon2D::segmentsIntersect(const Point2D& a1, const Point2D& b1,
                                  const Point2D& a2, const Point2D& b2) const {
    // Check if segments intersect using orientation tests
    auto orientation = [](const Point2D& p, const Point2D& q, const Point2D& r) -> double {
        return (q - p).cross(r - p);
    };
    
    double o1 = orientation(a1, b1, a2);
    double o2 = orientation(a1, b1, b2);
    double o3 = orientation(a2, b2, a1);
    double o4 = orientation(a2, b2, b1);
    if (o1 * o2 < 0 && o3 * o4 < 0) {
        return true;
    }
    
    // Special cases (collinear)
    auto onSegment = [](const Point2D& p, const Point2D& q, const Point2D& r) -> bool {
        return std::min(p.x, r.x) <= q.x && q.x <= std::max(p.x, r.x) &&
               std::min(p.y, r.y) <= q.y && q.y <= std::max(p.y, r.y);
    };
    
    if (std::abs(o1) < 1e-12 && onSegment(a1, a2, b1)) return true;
    if (std::abs(o2) < 1e-12 && onSegment(a1, b2, b1)) return true;
    if (std::abs(o3) < 1e-12 && onSegment(a2, a1, b2)) return true;
    if (std::abs(o4) < 1e-12 && onSegment(a2, b1, b2)) return true;
    
    return false;
}

double Polygon2D::perimeter() const {
    if (vertices.size() < 3) {
        return 0.0;
    }
    
    double perim = 0.0;
    for (size_t i = 0; i < vertices.size(); i++) {
        size_t j = (i + 1) % vertices.size();
        perim += vertices[i].distanceTo(vertices[j]);
    }
    
    return perim;
}

double pointToSegmentDistance(double px, double py, 
                              double x1, double y1, 
                              double x2, double y2) {
    double dx = x2 - x1;
    double dy = y2 - y1;
    double len2 = dx*dx + dy*dy;
    
    if (len2 < 1e-12) {
        return std::sqrt((px - x1)*(px - x1) + (py - y1)*(py - y1));
    }
    
    double t = ((px - x1)*dx + (py - y1)*dy) / len2;
    t = std::max(0.0, std::min(1.0, t));
    
    double cx = x1 + t * dx;
    double cy = y1 + t * dy;
    
    return std::sqrt((px - cx)*(px - cx) + (py - cy)*(py - cy));
}

bool pointInPolygon(double px, double py, const Polygon2D& poly) {
    bool inside = false;
    int n = poly.vertices.size();
    
    for (int i = 0, j = n - 1; i < n; j = i++) {
        double xi = poly.vertices[i].x, yi = poly.vertices[i].y;
        double xj = poly.vertices[j].x, yj = poly.vertices[j].y;
        
        bool intersect = ((yi > py) != (yj > py)) &&
                         (px < (xj - xi) * (py - yi) / (yj - yi) + xi);
        if (intersect) inside = !inside;
    }
    
    return inside;
}

double signedDistanceToPolygon(double px, double py, const Polygon2D& poly) {
    if (poly.vertices.empty()) return std::numeric_limits<double>::max();
    
    double minDist = std::numeric_limits<double>::max();
    int n = poly.vertices.size();
    
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        double dist = pointToSegmentDistance(px, py,
                                             poly.vertices[i].x, poly.vertices[i].y,
                                             poly.vertices[j].x, poly.vertices[j].y);
        minDist = std::min(minDist, dist);
    }
    
    bool inside = pointInPolygon(px, py, poly);
    return inside ? -minDist : minDist;
}