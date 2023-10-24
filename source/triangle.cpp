#include "triangle.h"

Point Triangle::bbox_centroid() const
{
    return (min(m_a, min(m_b, m_c)) + max(m_a, max(m_b, m_c))) / 2;
}

//From https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
bool Triangle::intersect(const Ray& ray, HitInfo& hit_info) const
{
    const float EPSILON = 0.0000001;
   Vector edge1, edge2, h, s, q;
   float a, f, u, v;
   edge1 = m_b - m_a;
   edge2 = m_c - m_a;

   h = cross(ray.direction, edge2);
   a = dot(edge1, h);

   if (a > -EPSILON && a < EPSILON)
       return false;    // This ray is parallel to this triangle.

   f = 1.0 / a;
   s = ray.origin - m_a;
   u = f * dot(s, h);

   if (u < 0.0 || u > 1.0)
       return false;

   q = cross(s, edge1);
   v = f * dot(ray.direction, q);

   if (v < 0.0 || u + v > 1.0)
       return false;

   // At this stage we can compute t to find out where the intersection point is on the line.
   float t = f * dot(edge2, q);

   if (t > EPSILON) // ray intersection
   {
       hit_info.inter_point = ray.origin + ray.direction * t;
       hit_info.normal_at_inter = normalize(cross(edge1, edge2));

       hit_info.t = t;

       hit_info.u = u;
       hit_info.v = v;
       return true;
   }
   else // This means that there is a line intersection but not a ray intersection.
       return false;
}

Point& Triangle::operator[] (int i)
{
    return *((&m_a) + i);
}

const Point& Triangle::operator[] (int i) const
{
    return *((&m_a) + i);
}
