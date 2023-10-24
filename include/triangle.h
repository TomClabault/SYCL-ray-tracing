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

    SYCL_EXTERNAL Point bbox_centroid() const;
    SYCL_EXTERNAL bool intersect(const Ray &ray, HitInfo& hit_info) const;

    SYCL_EXTERNAL Point& operator[] (int index);
    SYCL_EXTERNAL const Point& operator[] (int index) const;

    Point m_a, m_b, m_c;
};

#endif
