#include "intersections.h"

Eigen::Vector4d Point::pluC()
{
    Eigen::Vector4d coords;
    coords << r,1;
    return coords;
}

Eigen::Vector3d Point::operator-(const Point &p2)
{
    return (*this).r - p2.r;

}

Eigen::Vector3d Line::l()
{
    return r[1]-r[0];
}

Eigen::Vector3d Line::m()
{
    return r[0].cross(r[1]);
}

Vector6d Line::lm()
{
    Vector6d coords;
    coords << l(), m();
    return coords;
}

bool ifCoplanar(Line& l1, Line& l2)
{
    double shouldBeZero;
    double eps;
    eps = 10 * std::numeric_limits<double>::epsilon();
    shouldBeZero = l1.l().dot(l2.m())+l1.m().dot(l2.l());
    if (shouldBeZero <= eps) {
        return true;
    } else
        return false;
}

bool ifParallel(Line& l1, Line& l2)
{
    double shouldBeZero;
    double eps;
    eps = 10 * std::numeric_limits<double>::epsilon();
    shouldBeZero = l1.l().cross(l2.l()).norm();
    if (shouldBeZero <= eps) {
        return true;
    } else
        return false;
}

Eigen::Vector3d intersectionPoint(Line& l1, Line& l2)
{
    Eigen::Vector3d coords;
    if (ifCoplanar(l1, l2) & !ifParallel(l1, l2)) {
        Eigen::MatrixXd I(3,3);
        I.setIdentity();
        coords = (l1.l().dot(l2.m())*I + l1.l()*(l2.m().transpose()) - 
                l2.l()*(l1.m().transpose())) * l1.l().cross(l2.l()) 
                    / (l1.l().cross(l2.l())).squaredNorm();
    } else
        coords.setConstant(std::numeric_limits<double>::quiet_NaN());
    return coords;
}

Eigen::Vector3d Plane::n()
{
    return (r[1]-r[0]).cross((r[2]-r[0]));
}

Eigen::Vector3d Plane::n0()
{
    return n()/(n().norm());
}

double Plane::d0()
{
    return n0().dot(r[0]);
}

Eigen::Vector3d intersectionPointTUV(Line& line, Plane& plane)
{
    Eigen::Matrix<double, 3, 3> A;
    Eigen::Vector3d tuv;
    A << line.r[0] - line.r[1],
      plane.r[2] - plane.r[0],
      plane.r[1] - plane.r[0];
    tuv = A.inverse() * (line.r[0]-plane.r[0]);
    return tuv;
}

Eigen::Vector3d intersectionPoint(Line& line, Plane& plane)
{
    Eigen::Vector3d tuv, xyz;
    tuv = intersectionPointTUV(line, plane);
    xyz = line.r[0] + (line.r[1] - line.r[0]) * tuv(0);
    return xyz;
}

Eigen::Vector3d intersectionPoint(LineSegment& seg1, LineSegment& seg2)
{
    Line line1(seg1.r[0], seg1.r[1]), line2(seg2.r[0], seg2.r[1]);
    Eigen::Vector3d point;
    point = intersectionPoint(line1, line2);
    /*if (!isnan(point(0))) {
        if (!((point.array() >= seg1.r[0].cwiseMin(seg1.r[1]).array()).all() & 
            (point.array() <= seg1.r[0].cwiseMax(seg1.r[1]).array()).all() &
            (point.array() >= seg2.r[0].cwiseMin(seg2.r[1]).array()).all() & 
            (point.array() <= seg2.r[0].cwiseMax(seg2.r[1]).array()).all())) {
            point.setConstant(std::numeric_limits<double>::quiet_NaN());
        }
    }*/
    if (!(ifBelongsToLineSegment(point, seg1) & 
                ifBelongsToLineSegment (point, seg2))) {
        point.setConstant(std::numeric_limits<double>::quiet_NaN());
    } 
    return point;
}

Eigen::Vector3d intersectionPoint(LineSegment& segment, Triangle& tri)
{
    Eigen::Vector3d xyz, tuv;
    Plane plane(tri.r[0], tri.r[1], tri.r[2]);
    Line line(segment.r[0], segment.r[1]);
    tuv = intersectionPointTUV(line, plane);
    if ((tuv(0) >= 0) & (tuv(0) <= 1) &
            (tuv(1)>=0) & (tuv(2) >=0) & (tuv(1) + tuv(2) <=1)) {
        //    (tuv(1) + tuv(2) >=0) & (tuv(1) + tuv(2) <=1)) {
        xyz = line.r[0] + (line.r[1] - line.r[0]) * tuv(0);
        //xyz = intersectionPoint(line, plane);
    } else 
        xyz.setConstant(std::numeric_limits<double>::quiet_NaN());
    return xyz;
}

Eigen::Vector3d intersectionPoint(LineSegment& lSeg, Convex& cnvx)
{
    Triangle triangle;
    Eigen::Vector3d point;

    for (unsigned int i = 1 ; i < cnvx.vertices.size() - 1; i++) {
        triangle = Triangle(cnvx.vertices.at(0), cnvx.vertices.at(i),
                cnvx.vertices.at(i+1));
        point = intersectionPoint(lSeg, triangle);
        if (!std::isnan(point(0))) {break;}
    }
    return point;
}

std::vector<Eigen::Vector3d> intersectionPointN(LineSegment& lSeg, Convex& cnvx)
{
    /*There may be several crossings of such object. Function returns first found point.
     Errors may occur.*/
    Triangle triangle;
    Eigen::Vector3d point;
    std::vector<Eigen::Vector3d> rn;

    for (unsigned int i = 1 ; i < cnvx.vertices.size() - 1; i++) {
        triangle = Triangle(cnvx.vertices.at(0), cnvx.vertices.at(i),
                cnvx.vertices.at(i+1));
        point = intersectionPoint(lSeg, triangle);
        if (!std::isnan(point(0))) {break;}
    }
    rn.emplace_back(point);
    rn.emplace_back(triangle.n0());
    return rn;
}


bool ifBelongsToLine(Eigen::Vector3d& point, Line& line)
{
    double shouldBeZero;
    double eps; 
    eps = 10 * std::numeric_limits<double>::epsilon();
    shouldBeZero = (line.l().cross(point - line.r[0])).norm();
    if (shouldBeZero <= eps) {
        return true;
    } else
        return false;
}

bool ifBelongsToLineSegment(Eigen::Vector3d& point, LineSegment& lSeg)
{
    Line line(lSeg.r[0], lSeg.r[1]);
    if (ifBelongsToLine(point, line) & 
            (point.array() >= lSeg.r[0].cwiseMin(lSeg.r[1]).array()).all() & 
            (point.array() <= lSeg.r[0].cwiseMax(lSeg.r[1]).array()).all()) {
       return true;
    } else
       return false; 
}
