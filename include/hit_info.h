#ifndef HIT_INFO_H
#define HIT_INFO_H

#include "vec.h"

struct HitInfo
{
    Point inter_point;
    Vector normal_at_intersection;

    float t = -1.0f; //Distance along ray
    float u = -1, v = -1; //Barycentric coordinates
};

inline const sycl::stream operator << (const sycl::stream & os, const HitInfo & hit_info)
{
    os << "[" << "inter p: " << hit_info.inter_point << ", " << "normal: " << hit_info.normal_at_intersection << ", "  << ", t: " << hit_info.t << "]";

    return os;
}

#endif
