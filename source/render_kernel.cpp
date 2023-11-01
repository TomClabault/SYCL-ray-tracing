#include "render_kernel.h"

#include "triangle.h"
#include "vec.h"

Ray RenderKernel::get_camera_ray(float x, float y) const
{
    float x_ndc_space = x / m_width * 2 - 1;
    x_ndc_space *= (float)m_width / m_height; //Aspect ratio
    float y_ndc_space = y / m_height * 2 - 1;


    Point ray_origin_view_space(0.0f, 0.0f, 0.0f);
    Point ray_origin = m_camera.view_matrix(ray_origin_view_space);

    Point ray_point_direction_ndc_space = Point(x_ndc_space, y_ndc_space, m_camera.fov_dist);
    Point ray_point_direction_world_space = m_camera.view_matrix(ray_point_direction_ndc_space);

    Vector ray_direction = normalize(ray_point_direction_world_space - ray_origin);
    Ray ray(ray_origin, ray_direction);

    return ray;
}

void RenderKernel::ray_trace_pixel(int x, int y) const
{
    Color debug_final_color;
    Color final_color = Color(0.0f, 0.0f, 0.0f);
    for (int sample = 0; sample < SAMPLES_PER_KERNEL; sample++)
    {
        //Jittered around the center
        float x_jittered = (x + 0.5f);
        float y_jittered = (y + 0.5f);

        //TODO area sampling triangles
        Ray ray = get_camera_ray(x_jittered, y_jittered);

        Color throughput = Color(1.0f, 1.0f, 1.0f);
        Color sample_color = Color(0.0f, 0.0f, 0.0f);
        RayState next_ray_state = RayState::BOUNCE;

        for (int bounce = 0; bounce < MAX_BOUNCES; bounce++)
        {
            if (next_ray_state == BOUNCE)
            {
                HitInfo closest_hit_info;
                bool intersection_found = intersect_scene(ray, &closest_hit_info);

                if (intersection_found)
                {
                    sample_color = Color(abs(closest_hit_info.normal_at_intersection));

                    Ray shadow_ray(Point(0, 0, 0), Vector(1, 0, 0));

                    bool in_shadow = false;
                    in_shadow = evaluate_shadow_ray(shadow_ray, 2.0f);

                    if (in_shadow)
                        sample_color = Color();

                    Point new_ray_origin = closest_hit_info.inter_point + closest_hit_info.normal_at_intersection * 1.0e-4f;
                    Vector reflect_direction = normalize(ray.direction + 2 * dot(ray.direction, closest_hit_info.normal_at_intersection) * closest_hit_info.normal_at_intersection);
                    ray = Ray(new_ray_origin, reflect_direction);
                    next_ray_state = RayState::BOUNCE;
                }
                else
                    next_ray_state = RayState::MISSED;
            }
            else if (next_ray_state == MISSED)
                break;
            else if (next_ray_state == TERMINATED)
                break;
        }

        final_color += sample_color;
    }

    final_color /= SAMPLES_PER_KERNEL;
    final_color.a = 0.0f;
    m_frame_buffer_access[y * m_width + x] += final_color;
}

bool RenderKernel::intersect_scene(const Ray ray, HitInfo* closest_hit_info) const
{
    float closest_intersection_distance = -1.0f;
    bool inter_found = false;

    for (int i = 0; i < m_triangle_buffer_access.size(); i++)
    {
        const Triangle& triangle = m_triangle_buffer_access[i];

        HitInfo hit_info;
        if (triangle.intersect(ray, &hit_info))
        {
            if (hit_info.t < closest_intersection_distance || !inter_found)
            {
                closest_intersection_distance = hit_info.t;
                *closest_hit_info = hit_info;

                inter_found = true;
            }
        }
    }

    return inter_found;
}


bool RenderKernel::evaluate_shadow_ray(Ray ray, float t_max) const
{
    HitInfo local_hit_info;
    bool inter_found = intersect_scene(ray, &local_hit_info);

    if (inter_found)
    {
        if (local_hit_info.t + 1.0e-4f < t_max)
            //There is something in between the light and the origin of the ray
            return true;
        else
            return false;
    }

    return false;
}
