#ifndef SPHERE_H
#define SPHERE_H

#include "hit_info.h"
#include "ray.h"

struct Sphere
{
    Sphere(Point center, float radius, int material_index);

    SYCL_EXTERNAL bool intersect(const Ray& ray, HitInfo& hit_info) const;

    Point center;
    float radius;

    int material_index;
};

#endif
