#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "hit_info.h"
#include "ray.h"
#include "vec.h"

#include <sycl/sycl.hpp>

struct Triangle
{
    Triangle() {}
    Triangle(const Point& a, const Point& b, const Point& c) : m_a(a), m_b(b), m_c(c) {}

    SYCL_EXTERNAL bool intersect(Ray& ray, HitInfo& hit_info) const;

    Point m_a, m_b, m_c;
};

#endif
