//
//  eigenheaders.cpp
//  MultiTracker
//
//  Created by Fang Da on 10/27/14.
//
//

#include "eigenheaders.h"

Vec3i vc(const LosTopos::Vec3st & t)
{
    return Vec3i(t[0], t[1], t[2]);
}

Vec2i vc(const LosTopos::Vec2i& l)
{
    return Vec2i(l[0], l[1]);
}

Vec3d vc(const LosTopos::Vec3d & v)
{
    return Vec3d(v[0], v[1], v[2]);
}

LosTopos::Vec3d vc(const Vec3d & v)
{
    return LosTopos::Vec3d(v[0], v[1], v[2]);
}

