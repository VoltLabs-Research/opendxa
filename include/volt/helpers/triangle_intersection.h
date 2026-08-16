#pragma once

#include <volt/core/volt.h>

#include <array>
#include <cmath>

namespace Volt::TriangleIntersection{

using Triangle3 = std::array<Point3, 3>;

namespace detail{

inline double orient2d(const Point2& a, const Point2& b, const Point2& c) noexcept{
    return (a.x() - c.x()) * (b.y() - c.y()) - (a.y() - c.y()) * (b.x() - c.x());
}

inline bool testVertex2d(
    const Point2& p1, const Point2& q1, const Point2& r1,
    const Point2& p2, const Point2& q2, const Point2& r2
) noexcept{
    if(orient2d(r2, p2, q1) >= 0.0){
        if(orient2d(r2, q2, q1) <= 0.0){
            if(orient2d(p1, p2, q1) > 0.0){
                return orient2d(p1, q2, q1) <= 0.0;
            }
            if(orient2d(p1, p2, r1) < 0.0) return false;
            return orient2d(q1, r1, p2) >= 0.0;
        }
        if(orient2d(p1, q2, q1) > 0.0) return false;
        if(orient2d(r2, q2, r1) > 0.0) return false;
        return orient2d(q1, r1, q2) >= 0.0;
    }
    if(orient2d(r2, p2, r1) < 0.0) return false;
    if(orient2d(q1, r1, r2) >= 0.0){
        return orient2d(p1, p2, r1) >= 0.0;
    }
    if(orient2d(q1, r1, q2) < 0.0) return false;
    return orient2d(r2, r1, q2) >= 0.0;
}

inline bool testEdge2d(
    const Point2& p1, const Point2& q1, const Point2& r1,
    const Point2& p2, const Point2& q2, const Point2& r2
) noexcept{
    if(orient2d(r2, p2, q1) >= 0.0){
        if(orient2d(p1, p2, q1) >= 0.0){
            return orient2d(p1, q1, r2) >= 0.0;
        }
        if(orient2d(q1, r1, p2) < 0.0) return false;
        return orient2d(r1, p1, p2) >= 0.0;
    }
    if(orient2d(r2, p2, r1) < 0.0) return false;
    if(orient2d(p1, p2, r1) < 0.0) return false;
    if(orient2d(p1, r1, r2) >= 0.0) return true;
    return orient2d(q1, r1, r2) >= 0.0;
}

inline bool ccwOverlap2d(
    const Point2& p1, const Point2& q1, const Point2& r1,
    const Point2& p2, const Point2& q2, const Point2& r2
) noexcept{
    if(orient2d(p2, q2, p1) >= 0.0){
        if(orient2d(q2, r2, p1) >= 0.0){
            if(orient2d(r2, p2, p1) >= 0.0) return true;
            return testEdge2d(p1, q1, r1, p2, q2, r2);
        }
        if(orient2d(r2, p2, p1) >= 0.0){
            return testEdge2d(p1, q1, r1, r2, p2, q2);
        }
        return testVertex2d(p1, q1, r1, p2, q2, r2);
    }
    if(orient2d(q2, r2, p1) >= 0.0){
        if(orient2d(r2, p2, p1) >= 0.0){
            return testEdge2d(p1, q1, r1, q2, r2, p2);
        }
        return testVertex2d(p1, q1, r1, q2, r2, p2);
    }
    return testVertex2d(p1, q1, r1, r2, p2, q2);
}

inline bool overlap2d(
    const Point2& p1, const Point2& q1, const Point2& r1,
    const Point2& p2, const Point2& q2, const Point2& r2
) noexcept{
    const bool t1IsClockwise = orient2d(p1, q1, r1) < 0.0;
    const bool t2IsClockwise = orient2d(p2, q2, r2) < 0.0;

    if(t1IsClockwise){
        return t2IsClockwise
            ? ccwOverlap2d(p1, r1, q1, p2, r2, q2)
            : ccwOverlap2d(p1, r1, q1, p2, q2, r2);
    }
    return t2IsClockwise
        ? ccwOverlap2d(p1, q1, r1, p2, r2, q2)
        : ccwOverlap2d(p1, q1, r1, p2, q2, r2);
}

inline bool coplanarOverlap3d(
    const Point3& p1, const Point3& q1, const Point3& r1,
    const Point3& p2, const Point3& q2, const Point3& r2,
    const Vector3& normal1
) noexcept{
    const double nx = std::abs(normal1.x());
    const double ny = std::abs(normal1.y());
    const double nz = std::abs(normal1.z());

    if(nx > nz && nx >= ny){
        return overlap2d(
            {q1.z(), q1.y()}, {p1.z(), p1.y()}, {r1.z(), r1.y()},
            {q2.z(), q2.y()}, {p2.z(), p2.y()}, {r2.z(), r2.y()}
        );
    }
    if(ny > nz && ny >= nx){
        return overlap2d(
            {q1.x(), q1.z()}, {p1.x(), p1.z()}, {r1.x(), r1.z()},
            {q2.x(), q2.z()}, {p2.x(), p2.z()}, {r2.x(), r2.z()}
        );
    }
    return overlap2d(
        {p1.x(), p1.y()}, {q1.x(), q1.y()}, {r1.x(), r1.y()},
        {p2.x(), p2.y()}, {q2.x(), q2.y()}, {r2.x(), r2.y()}
    );
}

inline bool checkMinMax(
    const Point3& p1, const Point3& q1, const Point3& r1,
    const Point3& p2, const Point3& q2, const Point3& r2
) noexcept{
    if((q2 - q1).dot((p2 - q1).cross(p1 - q1)) > 0.0) return false;
    return (r2 - p1).dot((p2 - p1).cross(r1 - p1)) <= 0.0;
}

inline bool overlap3d(
    const Point3& p1, const Point3& q1, const Point3& r1,
    const Point3& p2, const Point3& q2, const Point3& r2,
    double dp2, double dq2, double dr2,
    const Vector3& normal1
) noexcept{
    if(dp2 > 0.0){
        if(dq2 > 0.0) return checkMinMax(p1, r1, q1, r2, p2, q2);
        if(dr2 > 0.0) return checkMinMax(p1, r1, q1, q2, r2, p2);
        return checkMinMax(p1, q1, r1, p2, q2, r2);
    }
    if(dp2 < 0.0){
        if(dq2 < 0.0) return checkMinMax(p1, q1, r1, r2, p2, q2);
        if(dr2 < 0.0) return checkMinMax(p1, q1, r1, q2, r2, p2);
        return checkMinMax(p1, r1, q1, p2, q2, r2);
    }
    if(dq2 < 0.0){
        if(dr2 >= 0.0) return checkMinMax(p1, r1, q1, q2, r2, p2);
        return checkMinMax(p1, q1, r1, p2, q2, r2);
    }
    if(dq2 > 0.0){
        if(dr2 > 0.0) return checkMinMax(p1, r1, q1, p2, q2, r2);
        return checkMinMax(p1, q1, r1, q2, r2, p2);
    }
    if(dr2 > 0.0) return checkMinMax(p1, q1, r1, r2, p2, q2);
    if(dr2 < 0.0) return checkMinMax(p1, r1, q1, r2, p2, q2);
    return coplanarOverlap3d(p1, q1, r1, p2, q2, r2, normal1);
}

}

[[nodiscard]] inline bool overlap(const Triangle3& t1, const Triangle3& t2) noexcept{
    const Point3& p1 = t1[0];
    const Point3& q1 = t1[1];
    const Point3& r1 = t1[2];
    const Point3& p2 = t2[0];
    const Point3& q2 = t2[1];
    const Point3& r2 = t2[2];

    const Vector3 normal2 = (p2 - r2).cross(q2 - r2);
    const double dp1 = (p1 - r2).dot(normal2);
    const double dq1 = (q1 - r2).dot(normal2);
    const double dr1 = (r1 - r2).dot(normal2);
    if(dp1 * dq1 > 0.0 && dp1 * dr1 > 0.0) return false;

    const Vector3 normal1 = (q1 - p1).cross(r1 - p1);
    const double dp2 = (p2 - r1).dot(normal1);
    const double dq2 = (q2 - r1).dot(normal1);
    const double dr2 = (r2 - r1).dot(normal1);
    if(dp2 * dq2 > 0.0 && dp2 * dr2 > 0.0) return false;

    if(dp1 > 0.0){
        if(dq1 > 0.0) return detail::overlap3d(r1, p1, q1, p2, r2, q2, dp2, dr2, dq2, normal1);
        if(dr1 > 0.0) return detail::overlap3d(q1, r1, p1, p2, r2, q2, dp2, dr2, dq2, normal1);
        return detail::overlap3d(p1, q1, r1, p2, q2, r2, dp2, dq2, dr2, normal1);
    }
    if(dp1 < 0.0){
        if(dq1 < 0.0) return detail::overlap3d(r1, p1, q1, p2, q2, r2, dp2, dq2, dr2, normal1);
        if(dr1 < 0.0) return detail::overlap3d(q1, r1, p1, p2, q2, r2, dp2, dq2, dr2, normal1);
        return detail::overlap3d(p1, q1, r1, p2, r2, q2, dp2, dr2, dq2, normal1);
    }
    if(dq1 < 0.0){
        if(dr1 >= 0.0) return detail::overlap3d(q1, r1, p1, p2, r2, q2, dp2, dr2, dq2, normal1);
        return detail::overlap3d(p1, q1, r1, p2, q2, r2, dp2, dq2, dr2, normal1);
    }
    if(dq1 > 0.0){
        if(dr1 > 0.0) return detail::overlap3d(p1, q1, r1, p2, r2, q2, dp2, dr2, dq2, normal1);
        return detail::overlap3d(q1, r1, p1, p2, q2, r2, dp2, dq2, dr2, normal1);
    }
    if(dr1 > 0.0) return detail::overlap3d(r1, p1, q1, p2, q2, r2, dp2, dq2, dr2, normal1);
    if(dr1 < 0.0) return detail::overlap3d(r1, p1, q1, p2, r2, q2, dp2, dr2, dq2, normal1);
    return detail::coplanarOverlap3d(p1, q1, r1, p2, q2, r2, normal1);
}

}
