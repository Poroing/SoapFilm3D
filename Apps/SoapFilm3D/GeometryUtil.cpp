#include "GeometryUtil.h"

double
angleAroundAxis(const Vec3d& v0,
                const Vec3d& v1,
                const Vec3d& a) // angle from v0 to v1 around axis a
{
    double asq = a.squaredNorm();
    assert(asq != 0);

    Vec3d u = v0 - v0.dot(a) * a / asq;
    assert(u.squaredNorm() != 0);
    u.normalize();

    Vec3d v = a.cross(u).normalized();

    return atan2(v1.dot(v), v1.dot(u));
}
