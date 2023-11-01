
#ifndef _MAT_H
#define _MAT_H

#include "vec.h"

struct Transform
{
    SYCL_EXTERNAL Transform(
        const float t00 = 1.0f, const float t01 = 0.0f, const float t02 = 0.0f, const float t03 = 0.0f,
        const float t10 = 0.0f, const float t11 = 1.0f, const float t12 = 0.0f, const float t13 = 0.0f,
        const float t20 = 0.0f, const float t21 = 0.0f, const float t22 = 1.0f, const float t23 = 0.0f,
        const float t30 = 0.0f, const float t31 = 0.0f, const float t32 = 0.0f, const float t33 = 1.0f)
    {
        m[0][0]= t00; m[0][1]= t01; m[0][2]= t02; m[0][3]= t03;
        m[1][0]= t10; m[1][1]= t11; m[1][2]= t12; m[1][3]= t13;
        m[2][0]= t20; m[2][1]= t21; m[2][2]= t22; m[2][3]= t23;
        m[3][0]= t30; m[3][1]= t31; m[3][2]= t32; m[3][3]= t33;
    }

    SYCL_EXTERNAL Transform(const Vector& x, const Vector& y, const Vector& z, const Vector& w)
    {
        m[0][0] = x.x;	m[0][1] = y.x;	m[0][2] = z.x;	m[0][3] = w.x;
        m[1][0] = x.y;	m[1][1] = y.y;	m[1][2] = z.y;	m[1][3] = w.y;
        m[2][0] = x.z;	m[2][1] = y.z;	m[2][2] = z.z;	m[2][3] = w.z;
        m[3][0] = 0;	m[3][1] = 0;	m[3][2] = 0;	m[3][3] = 1;
    }

    SYCL_EXTERNAL Vector operator[] (const unsigned c) const;

    SYCL_EXTERNAL Point operator() (const Point& p) const;
    SYCL_EXTERNAL Vector operator() (const Vector& v) const;

    SYCL_EXTERNAL const float* data() const { return &m[0][0]; }

    float m[4][4];
};

SYCL_EXTERNAL Transform Identity();

SYCL_EXTERNAL Transform Translation(const Vector& v);
SYCL_EXTERNAL Transform Translation(const float x, const float y, const float z);

SYCL_EXTERNAL Transform compose_transform(const Transform& a, const Transform& b);
SYCL_EXTERNAL Transform operator* (const Transform& a, const Transform& b);

#endif
