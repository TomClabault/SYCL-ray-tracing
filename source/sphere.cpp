#include "sphere.h"

Sphere::Sphere(Point center, float radius, int material_index) : center(center), radius(radius), material_index(material_index) { }

bool Sphere::intersect(const Ray &ray, HitInfo& hit_info) const
{
    Vector L = ray.origin - center;

    //dot(ray._direction, ray._direction) = 1 because direction is normalized
    constexpr float a = 1;
    float b = 2 * dot(ray.direction, L);
    float c = dot(L, L) - radius * radius;

    float delta = b * b - 4 * a * c;
    if (delta < 0)
        return false;
    else
    {
        constexpr float a2 = 2 * a;

        if (delta == 0.0)
            hit_info.t = -b / a2;
        else
        {
            float sqrt_delta = std::sqrt(delta);

            float t1 = (-b - sqrt_delta) / a2;
            float t2 = (-b + sqrt_delta) / a2;

            if (t1 < t2)
            {
                hit_info.t = t1;
                if (hit_info.t < 0)
                    hit_info.t = t2;
            }
        }

        if (hit_info.t < 0)
            return false;

        hit_info.normal_at_intersection = normalize((ray.origin + ray.direction * hit_info.t) - center);
        hit_info.material_index = material_index;

        return true;
    }
}
