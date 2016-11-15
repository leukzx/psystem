#include <iostream>
#include "intersections.h"

int main(int argc, char **argv)
{
    Eigen::Vector3d a(1,0,0);
    Eigen::Vector3d b(0,1,0);
    //Eigen::Vector3d c(3,1,0);
    //Eigen::Vector3d d(1,2,0);
    Eigen::Vector3d c(0,0,1);
    Eigen::Vector3d d(1,1,1);
    Eigen::Vector3d e(-2,-2,-2);
    Eigen::Vector3d f(3,2,1);

    Line L1(a,e);
    Line L2(d,e);

    Plane P(c,a,b);

    std::cout << L1.lm().transpose() << std::endl;
    std::cout << L2.lm().transpose() << std::endl;
    std::cout << "Coplanar: " << ifCoplanar(L1, L2) << std::endl;
    std::cout << "Parallel: " << ifParallel(L1, L2) << std::endl;
    std::cout << "L1 & L2 intersect at: " << intersectionPoint(L1, L2).transpose() << std::endl;
    std::cout << "L2 & P intersect at:\n " << intersectionPoint(L2, P).transpose() << std::endl;
    std::cout << P.d0() << std::endl;

    std::cout << P.n0() << std::endl;
    std::cout << P.n() << std::endl;

    LineSegment lSeg(d, e);
    Triangle Tri(a,b,c);

    std::cout << "Triangle - LineSegment intersection point: "
              << intersectionPoint(lSeg, Tri).transpose() << std::endl;
    std::cout << "Tiangle\'s normal vector: " << Tri.n0().transpose() << std::endl;


    std::vector<Eigen::Vector3d> convexPoints;
    convexPoints.emplace_back(a);
    convexPoints.emplace_back(b);
    convexPoints.emplace_back(c);
    Convex convex(convexPoints);
    std::cout << "Convex - LineSegment intersection point: "
              << intersectionPoint(lSeg, convex).transpose() << std::endl;
    
    LineSegment lSeg2(a, f);
    std::cout << "lSeg - lSeg2 intersection: " << intersectionPoint(lSeg, lSeg2).transpose() << std::endl;
    
    Eigen::Vector3d g(-1, 0, 0);
    Eigen::Vector3d h(-2, 0, 0);
    Eigen::Vector3d i(1, 0, 0);
    Eigen::Vector3d j(0, -2, 0);
    Eigen::Vector3d k(0, 1, 0);
    Line line1(h, i);
    Line line2(j, k);
    LineSegment lineSegm1(h, i);
    LineSegment lineSegm2(j, k);
    std::cout << "g belongs to line1: " << ifBelongsToLine(g, line1) << std::endl;
    std::cout << "g belongs to lineSegm1: " << ifBelongsToLineSegment(g, lineSegm1) << std::endl;
    std::cout << "line1 & line2 intersection point is: " << intersectionPoint(line1, line2).transpose() << std::endl;
    std::cout << "lineSegm1 & lineSegm2 intersection point is: " << intersectionPoint(lineSegm1, lineSegm2).transpose() << std::endl;
}
