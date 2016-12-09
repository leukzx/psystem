#include <iostream>
#include <cmath>
//#include <eigen3/Eigen/Dense>
#include "/usr/include/eigen3/Eigen/Dense"
#include <fstream>      // Stream class to both read and write from/to files.
#include <string>       // std::string, std::to_string
#include <iomanip>      // std::setprecision
#include <cfloat>       // DBL_DIG (for output precision)
#include <vector>

typedef Eigen::Matrix<double, 6, 1> Vector6d;

class Line;
class Plane;
class LineSegment;
class Triangle;
class Convex;

class Point {
 public:
    Eigen::Vector3d r;
    
    Eigen::Vector4d pluC(); //Plucker coordinates

    Eigen::Vector3d operator-(const Point&);

    Point() {r << 0,0,0;}
    Point(Eigen::Vector3d R): r(R) {}
    Point(double x, double y, double z) {r << x,y,z;}
};

class Line {
 public:

    std::vector<Eigen::Vector3d> r;

    Eigen::Vector3d l(); //Direction vector
    Eigen::Vector3d m(); //Moment vector about the origin
    Vector6d lm(); //Plucker coordinates for line
    
    Line() {};
    Line(Eigen::Vector3d &p0, Eigen::Vector3d &p1): r({p0, p1}) {}
};

class Plane {
  public:
    std::vector<Eigen::Vector3d> r;

    Eigen::Vector3d n(); //Normal vector
    Eigen::Vector3d n0(); //Normalized normal vector
    double d0(); //Distance from the origin
    Plane() {};
    Plane(Eigen::Vector3d &p0, Eigen::Vector3d &p1, Eigen::Vector3d &p2):
        r({p0, p1, p2}) {}

};

class LineSegment: public Line {
  public:
    LineSegment() {};  
    LineSegment(Eigen::Vector3d &p0, Eigen::Vector3d &p1): Line(p0, p1) {}
};

class Triangle: public Plane {
 public:
    Triangle() {};
    Triangle(Eigen::Vector3d &p0, Eigen::Vector3d &p1, Eigen::Vector3d &p2):
        Plane(p0, p1, p2) {}
};

class Convex {
  public:
    std::vector<Eigen::Vector3d> vertices;
    Convex() {};
    Convex(std::vector<Eigen::Vector3d>& v): vertices(v) {}
};

//Returns true if point belongs to line
bool ifBelongsToLine(Eigen::Vector3d& point, Line& line); 
//Returns true if point belongs to line segment
bool ifBelongsToLineSegment(Eigen::Vector3d& point, LineSegment& lSegment); 
bool ifCoplanar(Line& l1, Line& l2);
bool ifParallel(Line& l1, Line& l2);
Eigen::Vector3d intersectionPoint(Line& l1, Line& l2);
Eigen::Vector3d intersectionPoint(Line& line, Plane& plane);
Eigen::Vector3d intersectionPoint(LineSegment& ls1, LineSegment& ls2);
Eigen::Vector3d intersectionPoint(LineSegment& segment, Triangle& tri);
Eigen::Vector3d intersectionPoint(LineSegment& lSeg, Convex& cnvx);
std::vector<Eigen::Vector3d> intersectionPointN(LineSegment& lSeg, Convex& cnvx); // Returns cross point and normal vector
