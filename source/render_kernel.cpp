#include "render_kernel.h"

#include "triangle.h"

void RenderKernel::operator()(const sycl::nd_item<2>& coordinates) const
{
    int x = coordinates.get_global_id(0);
    int y = coordinates.get_global_id(1);

    ray_trace_pixel(x, y);
}

void branchlessONB(const Vector& n, Vector& b1, Vector& b2)
{
    float sign = sycl::copysign(1.0f, n.z);
    const float a = -1.0f / (sign + n.z);
    const float b = n.x * n.y * a;
    b1 = Vector(1.0f + sign * n.x * n.x * a, sign * b, -sign * n.x);
    b2 = Vector(b, sign + n.y * n.y * a, -n.y);
}

Vector RenderKernel::rotate_vector_around_normal(const Vector& normal, const Vector& random_dir_local_space) const
{
    Vector tangent, bitangent;
    branchlessONB(normal, tangent, bitangent);

    //Transforming from the random_direction in its local space to the space around the normal
    //given in parameter (the space with the given normal as the Z up vector)
    return random_dir_local_space.x * tangent + random_dir_local_space.y * bitangent + random_dir_local_space.z * normal;
}

Vector RenderKernel::uniform_direction_around_normal(const Vector& normal, float& pdf, xorshift32_generator& random_number_generator) const
{
    float rand_1 = random_number_generator();
    float rand_2 = random_number_generator();

    float phi = 2.0f * M_PI * rand_1;
    float root = sycl::sqrt(1 - rand_2 * rand_2);

    pdf = 1.0f / (2.0f * M_PI);

    //Generating a random direction in a local space with Z as the Up vector
    Vector random_dir_local_space(sycl::cos(phi) * root, sycl::sin(phi) * root, rand_2);
    return rotate_vector_around_normal(normal, random_dir_local_space);
}

Vector RenderKernel::cosine_weighted_direction_around_normal(const Vector& normal, float& pdf, xorshift32_generator& random_number_generator) const
{
    float rand_1 = random_number_generator();
    float rand_2 = random_number_generator();

    float sqrt_rand_2 = sycl::sqrt(rand_2);
    float phi = 2.0f * M_PI * rand_1;
    float cos_theta = sqrt_rand_2;
    float sin_theta = sycl::sqrt(sycl::max(0.0f, 1.0f - cos_theta * cos_theta));

    pdf = sqrt_rand_2 / M_PI;

    //Generating a random direction in a local space with Z as the Up vector
    Vector random_dir_local_space = Vector(sycl::cos(phi) * sin_theta, sycl::sin(phi) * sin_theta, sqrt_rand_2);
    return rotate_vector_around_normal(normal, random_dir_local_space);
}

Ray RenderKernel::get_camera_ray(float x, float y) const
{
    float x_ndc_space = x / m_width * 2 - 1;
    x_ndc_space *= (float)m_width / m_height; //Aspect ratio
    float y_ndc_space = y / m_height * 2 - 1;


    Point ray_origin_view_space(0, 0, 0);
    Point ray_origin = m_camera.view_matrix(ray_origin_view_space);

    Point ray_point_direction_ndc_space = Point(x_ndc_space, y_ndc_space, m_camera.fov_dist);
    Point ray_point_direction_world_space = m_camera.view_matrix(ray_point_direction_ndc_space);

    Vector ray_direction = normalize(ray_point_direction_world_space - ray_origin);
    Ray ray(ray_origin, ray_direction);

    return ray;
}

void RenderKernel::ray_trace_pixel(int x, int y) const
{
    xorshift32_generator random_number_generator(x * y * SAMPLES_PER_KERNEL * (m_kernel_iteration + 1));

    Color final_color = Color(0.0f, 0.0f, 0.0f);
    for (int sample = 0; sample < SAMPLES_PER_KERNEL; sample++)
    {
        //Jittered around the center
        float x_jittered = (x + 0.5f) + random_number_generator() - 1.0f;
        float y_jittered = (y + 0.5f) + random_number_generator() - 1.0f;

        //TODO area sampling triangles
        Ray ray = get_camera_ray(x_jittered, y_jittered);

        Color throughput = Color(1.0f, 1.0f, 1.0f);
        Color sample_color = Color(0.0f, 0.0f, 0.0f);
        RayState next_ray_state = RayState::BOUNCE;
        float random_direction_pdf = 1.0f;
        for (int bounce = 0; bounce < MAX_BOUNCES; bounce++)
        {
            if (next_ray_state == BOUNCE)
            {
                HitInfo closest_hit_info;
                bool intersection_found = intersect_scene(ray, closest_hit_info);

                if (intersection_found)
                {
                    //Indirect lighting
                    int material_index = m_materials_indices_buffer[closest_hit_info.triangle_index];
                    SimpleMaterial mat = m_materials_buffer_access[material_index];

                    throughput *= mat.diffuse * sycl::max(0.0f, dot(-ray.direction, closest_hit_info.normal_at_inter));
                    if (bounce == 0)
                        sample_color += mat.emission;


                    // Direct lighting
                    float pdf;
                    LightSourceInformation light_source_info;
                    Point random_light_point = sample_random_point_on_lights(random_number_generator, pdf, light_source_info);
                    Point shadow_ray_origin = closest_hit_info.inter_point + closest_hit_info.normal_at_inter * 1.0e-4f;
                    Vector shadow_ray_direction = random_light_point - shadow_ray_origin;
                    float distance_to_light = length(shadow_ray_direction);
                    Vector shadow_ray_direction_normalized = normalize(shadow_ray_direction);

                    Ray shadow_ray(shadow_ray_origin, shadow_ray_direction_normalized);

                    bool in_shadow = evaluate_shadow_ray(shadow_ray, distance_to_light);

                    Color radiance = Color(0.0f, 0.0f, 0.0f);
                    if (!in_shadow)
                    {
                        const SimpleMaterial& emissive_triangle_material = m_materials_buffer_access[m_materials_indices_buffer[light_source_info.emissive_triangle_index]];

                        radiance = emissive_triangle_material.emission;
                        //Cosine angle on the illuminated surface
                        radiance *= sycl::max(dot(closest_hit_info.normal_at_inter, shadow_ray_direction_normalized), 0.0f);
                        //Cosine angle on the light surface
                        radiance *= sycl::max(dot(light_source_info.light_source_normal, -shadow_ray_direction_normalized), 0.0f);
                        //Falloff of the light intensity with the distance squared
                        radiance /= distance_to_light * distance_to_light;
                        //PDF: Probability of having chosen this point on this exact light source
                        radiance /= pdf;
                        //The illuminated surface is Lambertian
                        radiance /= M_PI;
                    }

                    Vector random_dir = uniform_direction_around_normal(closest_hit_info.normal_at_inter, random_direction_pdf, random_number_generator);
                    Point new_ray_origin = closest_hit_info.inter_point + closest_hit_info.normal_at_inter * 1.0e-4f;

                    ray = Ray(new_ray_origin, normalize(random_dir));
                    next_ray_state = RayState::BOUNCE;

                    sample_color += radiance * throughput;

                    //Cosine angle of the bounced ray
                    throughput /= random_direction_pdf;
                }
                else
                    next_ray_state = RayState::MISSED;
            }
            else if (next_ray_state == MISSED)
            {
                //Handle skysphere here
                break;
            }
            else if (next_ray_state == TERMINATED)
                break;
        }

        final_color += sample_color;
    }

    final_color /= SAMPLES_PER_KERNEL;
    final_color.a = 0.0f;
    m_frame_buffer_access[y * m_width + x] += final_color;

    if (m_kernel_iteration == RENDER_KERNEL_ITERATIONS - 1)
    {
        //Last iteration, computing the average
        m_frame_buffer_access[y * m_width + x] /= RENDER_KERNEL_ITERATIONS;

        const float gamma = 2.2;
        const float exposure = 1.0f;
        Color hdrColor = m_frame_buffer_access[y * m_width + x];

        //Exposure tone mapping
        Color tone_mapped = Color(1.0f, 1.0f, 1.0f) - exp(-hdrColor * exposure);
        // Gamma correction
        Color gamma_corrected = pow(tone_mapped, 1.0f / gamma);

        m_frame_buffer_access[y * m_width + x] = gamma_corrected;
    }
}

bool RenderKernel::intersect_scene(Ray& ray, HitInfo& closest_hit_info) const
{
    float closest_intersection_distance = -1;
    bool intersection_found = false;

    for (int i = 0; i < m_triangle_buffer_access.size(); i++)
    {
        const Triangle triangle = m_triangle_buffer_access[i];

        HitInfo hit_info;
        if (triangle.intersect(ray, hit_info))
        {
            if (hit_info.t < closest_intersection_distance || closest_intersection_distance == -1.0f)
            {
                hit_info.triangle_index = i;
                closest_intersection_distance = hit_info.t;
                closest_hit_info = hit_info;

                intersection_found = true;
            }
        }
    }

    return intersection_found;
}

Point RenderKernel::sample_random_point_on_lights(xorshift32_generator& random_number_generator, float& pdf, LightSourceInformation& light_info) const
{
    light_info.emissive_triangle_index = random_number_generator() * m_emissive_triangle_indices_buffer.size();
    light_info.emissive_triangle_index = m_emissive_triangle_indices_buffer[light_info.emissive_triangle_index];
    Triangle random_emissive_triangle = m_triangle_buffer_access[light_info.emissive_triangle_index];

    float rand_1 = random_number_generator();
    float rand_2 = random_number_generator();

    float sqrt_r1 = sycl::sqrt(rand_1);
    float u = 1.0f - sqrt_r1;
    float v = (1.0f - rand_2) * sqrt_r1;

    Vector AB = random_emissive_triangle.m_b - random_emissive_triangle.m_a;
    Vector AC = random_emissive_triangle.m_c - random_emissive_triangle.m_a;

    Point random_point_on_triangle = random_emissive_triangle.m_a + AB * u + AC * v;

    Vector normal = cross(AB, AC);
    float length_normal = length(normal);
    light_info.light_source_normal = normal / length_normal; //Normalized
    float triangle_area = length_normal * 0.5f;
    float nb_emissive_triangles = m_emissive_triangle_indices_buffer.size();

    pdf = 1.0f / (nb_emissive_triangles * triangle_area);

    return random_point_on_triangle;
}

bool RenderKernel::evaluate_shadow_ray(Ray& ray, float t_max) const
{
    HitInfo hit_info;
    intersect_scene(ray, hit_info);
    if (hit_info.t + 1.0e-4f < t_max)
        //There is something in between the light and the origin of the ray
        return true;
    else
        return false;
}
