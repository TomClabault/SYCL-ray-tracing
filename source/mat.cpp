#define _USE_MATH_DEFINES

#include <cassert>
#include <cstdio>
#include <cmath>
#include <algorithm>

#include "mat.h"

Vector Transform::operator[] ( const unsigned c ) const
{
    assert(c < 4);
    return Vector(m[0][c], m[1][c], m[2][c]);
}


Point Transform::operator() ( const Point& p ) const
{
    float x= p.x;
    float y= p.y;
    float z= p.z;

    float xt= m[0][0] * x + m[0][1] * y + m[0][2] * z + m[0][3];        // dot(vec4(m[0]), vec4(p, 1))
    float yt= m[1][0] * x + m[1][1] * y + m[1][2] * z + m[1][3];        // dot(vec4(m[1]), vec4(p, 1))
    float zt= m[2][0] * x + m[2][1] * y + m[2][2] * z + m[2][3];        // dot(vec4(m[2]), vec4(p, 1))
    float wt= m[3][0] * x + m[3][1] * y + m[3][2] * z + m[3][3];        // dot(vec4(m[3]), vec4(p, 1))

    assert(wt != 0);
    float w= 1.f / wt;
    if(wt == 1.f)
        return Point(xt, yt, zt);
    else
        return Point(xt*w, yt*w, zt*w);
}

Vector Transform::operator() ( const Vector& v ) const
{
    float x= v.x;
    float y= v.y;
    float z= v.z;

    float xt= m[0][0] * x + m[0][1] * y + m[0][2] * z;                  // dot(vec4(m[0]), vec4(v, 0))
    float yt= m[1][0] * x + m[1][1] * y + m[1][2] * z;                  // dot(vec4(m[1]), vec4(v, 0))
    float zt= m[2][0] * x + m[2][1] * y + m[2][2] * z;                  // dot(vec4(m[2]), vec4(v, 0))
    // dot(vec4(m[3]), vec4(v, 0)) == dot(vec4(0, 0, 0, 1), vec4(v, 0)) == 0 by definition

    return Vector(xt, yt, zt);
}

Transform Identity( )
{
    return Transform();
}

Transform Translation( const Vector& v )
{
    return Transform(
        1, 0, 0, v.x,
        0, 1, 0, v.y,
        0, 0, 1, v.z,
        0, 0, 0, 1);
}

Transform Translation( const float x, const float y, const float z )
{
    return Translation( Vector(x, y, z) );
}

Transform compose_transform( const Transform& a, const Transform& b )
{
    Transform m;
    for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
            m.m[i][j]= a.m[i][0] * b.m[0][j] + a.m[i][1] * b.m[1][j] + a.m[i][2] * b.m[2][j] + a.m[i][3] * b.m[3][j];

    return m;
}

Transform operator* ( const Transform& a, const Transform& b )
{
    return compose_transform(a, b);
}
