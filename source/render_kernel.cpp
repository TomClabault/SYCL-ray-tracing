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
    if (x != m_width / 2 || y != m_height / 2)
        return;

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
        for (int bounce = 0; bounce < MAX_BOUNCES; bounce++)
        {
            if (next_ray_state == BOUNCE)
            {
                HitInfo closest_hit_info;
                bool intersection_found = intersect_scene(ray, closest_hit_info);

                if (intersection_found)
                {
                    int material_index = m_materials_indices_buffer[closest_hit_info.triangle_index];
                    SimpleMaterial material = m_materials_buffer_access[material_index];

                    // ------------------------------------- //
                    // ---------- Direct lighting ---------- //
                    // ------------------------------------- //
                    float light_sample_pdf;
                    LightSourceInformation light_source_info;
                    Point random_light_point = sample_random_point_on_lights(random_number_generator, light_sample_pdf, light_source_info);
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
                        radiance /= light_sample_pdf;
                        //BRDF of the illuminated surface
                        radiance *= cook_torrance_brdf(material, shadow_ray.direction, -ray.direction, closest_hit_info.normal_at_inter);
                    }

                    float random_bounce_direction_pdf;
                    //Vector random_bounce_direction = uniform_direction_around_normal(closest_hit_info.normal_at_inter, random_bounce_direction_pdf, random_number_generator);



                    // --------------------------------------- //
                    // ---------- Indirect lighting ---------- //
                    // --------------------------------------- //

                    //TODO dans quel sens on doit appliquer le throughput, pdf, ...
                    Vector random_bounce_direction;
                    Point new_ray_origin = closest_hit_info.inter_point + closest_hit_info.normal_at_inter * 1.0e-4f;
                    Color brdf = cook_torrance_brdf_importance_sample(material, -ray.direction, closest_hit_info.normal_at_inter, random_bounce_direction, random_bounce_direction_pdf, random_number_generator);
                    //Color brdf = cook_torrance_brdf(material, random_bounce_direction, -ray.direction, closest_hit_info.normal_at_inter);
                    throughput *= brdf * sycl::max(0.0f, dot(random_bounce_direction, closest_hit_info.normal_at_inter));

                    if (bounce == 0)
                        sample_color += material.emission;
                    sample_color += radiance * throughput;

                    throughput /= random_bounce_direction_pdf;

                    ray = Ray(new_ray_origin, random_bounce_direction);
                    next_ray_state = RayState::BOUNCE;
                }
                else
                    next_ray_state = RayState::MISSED;
            }
            else if (next_ray_state == MISSED)
            {
                float u, v;
                u = 0.5f + sycl::atan2(ray.direction.z, ray.direction.x) / (2.0f * (float)M_PI);
                v = 0.5f + sycl::asin(ray.direction.y) / (float)M_PI;

                sycl::int2 coords(v * m_skysphere_height, u * m_skysphere_width);
                sycl::float4 skysphere_color_float4 = m_skysphere.read(coords, m_skysphere_sampler);
                Color skysphere_color = Color(skysphere_color_float4.x(), skysphere_color_float4.y(), skysphere_color_float4.z());

                //sample_color += skysphere_color * throughput;

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
        const float exposure = 2.5f;
        Color hdrColor = m_frame_buffer_access[y * m_width + x];

        //Exposure tone mapping
        Color tone_mapped = Color(1.0f, 1.0f, 1.0f) - exp(-hdrColor * exposure);
        // Gamma correction
        Color gamma_corrected = pow(tone_mapped, 1.0f / gamma);

        m_frame_buffer_access[y * m_width + x] = gamma_corrected;
    }
}

Color RenderKernel::lambertian_brdf(const SimpleMaterial& material, const Vector& to_light_direction, const Vector& view_direction, const Vector& surface_normal) const
{
    return material.diffuse / M_PI;
}

Color fresnel_schlick(Color F0, float NoV)
{
    return F0 + (Color(1.0f) - F0) * sycl::pow((1.0f - NoV), 5.0f);
}

float GGX_normal_distribution(float alpha, float NoH)
{
    float alpha2 = alpha * alpha;
    float NoH2 = NoH * NoH;
    float b = (NoH2 * (alpha2 - 1.0f) + 1.0f);
    return alpha2 * (1.0f / M_PI) / (b * b);
}

float G1_schlick_ggx(float k, float dot_prod)
{
    return dot_prod / (dot_prod * (1.0f - k) + k);
}

float GGX_smith_masking_shadowing(float roughness_squared, float NoV, float NoL)
{
    float k = roughness_squared / 2.0f;

    return G1_schlick_ggx(k, NoL) * G1_schlick_ggx(k, NoV);
}

Color RenderKernel::cook_torrance_brdf(const SimpleMaterial& material, const Vector& to_light_direction, const Vector& view_direction, const Vector& surface_normal) const
{
    Color brdf_color = Color(0.0f, 0.0f, 0.0f);
    Color base_color = material.diffuse;

    Vector halfway_vector = normalize(view_direction + to_light_direction);

    float NoV = sycl::max(0.0f, dot(surface_normal, view_direction));
    float NoL = sycl::max(0.0f, dot(surface_normal, to_light_direction));
    float NoH = sycl::max(0.0f, dot(surface_normal, halfway_vector));
    float VoH = sycl::max(0.0f, dot(halfway_vector, view_direction));

    if (NoV > 0.0f && NoL > 0.0f && NoH > 0.0f)
    {
        float metalness = 0.0f;//texture2D(u_mesh_specular_texture, vs_texcoords).b;
        float roughness = 0.1f;//texture2D(u_mesh_specular_texture, vs_texcoords).g;

        float alpha = roughness * roughness;

        ////////// Cook Torrance BRDF //////////
        Color F;
        float D, G;

        //F0 = 0.04 for dielectrics, 1.0 for metals (approximation)
        Color F0 = Color(0.04f * (1.0f - metalness)) + metalness * base_color;

        //GGX Distribution function
        F = fresnel_schlick(F0, VoH);
        D = GGX_normal_distribution(alpha, NoH);
        G = GGX_smith_masking_shadowing(alpha, NoV, NoL);

        Color kD = Color(1.0f - metalness); //Metals do not have a diffuse part
        kD *= Color(1.0f) - F;//Only the transmitted light is diffused

        Color diffuse_part = kD * base_color / M_PI;
        Color specular_part = (F * D * G) / (4.0f * NoV * NoL);

        brdf_color = diffuse_part + specular_part;
    }

    return brdf_color;
}

Color RenderKernel::cook_torrance_brdf_importance_sample(const SimpleMaterial& material, const Vector& view_direction, const Vector& surface_normal, Vector& output_direction, float& pdf, xorshift32_generator& random_number_generator) const
{
    float metalness = 0.0f;
    float roughness = 0.1f;

    float alpha = roughness * roughness;

    float rand1 = random_number_generator();
    float rand2 = random_number_generator();

    float phi = 2.0f * M_PI * rand1;
    float theta = sycl::acos((1 - rand2) / (rand2 * (alpha * alpha - 1) + 1));
    float sin_theta = sycl::sin(theta);

    Vector microfacet_normal_local_space = normalize(Vector(sycl::cos(phi) * sin_theta, sycl::sin(phi) * sin_theta, sycl::cos(theta)));
    //m_out_stream << microfacet_normal_local_space << sycl::endl;
    Vector microfacet_normal = normalize(rotate_vector_around_normal(surface_normal, microfacet_normal_local_space));//TODO remove normalize ?
    if (dot(microfacet_normal, surface_normal) < 0)
        //The microfacet normal that we sampled was under the surface, it can happen
        return Color(0, 0, 0);
    Vector to_light_direction = normalize(2.0f * dot(microfacet_normal, view_direction) * microfacet_normal - view_direction);//TODO remove normalize ?
    output_direction = to_light_direction;
    Vector halfway_vector = microfacet_normal;
    //halfway_vector = normalize(surface_normal + view_direction);

    Color brdf_color = Color(0.0f, 0.0f, 0.0f);
    Color base_color = material.diffuse;

    float NoV = sycl::max(0.0f, dot(surface_normal, view_direction));
    float NoL = sycl::max(0.0f, dot(surface_normal, to_light_direction));
    float NoH = sycl::max(0.0f, dot(surface_normal, halfway_vector));
    float VoH = sycl::max(0.0f, dot(halfway_vector, view_direction));

    if (NoV > 0.0f && NoL > 0.0f && NoH > 0.0f)
    {
        /////////// Cook Torrance BRDF //////////
        Color F;
        float D, G;

        //F0 = 0.04 for dielectrics, 1.0 for metals (approximation)
        Color F0 = Color(0.04f * (1.0f - metalness)) + metalness * base_color;

        //GGX Distribution function
        F = fresnel_schlick(F0, VoH);
        D = GGX_normal_distribution(alpha, NoH);
        m_out_stream << D << sycl::endl;
        G = GGX_smith_masking_shadowing(alpha, NoV, NoL);

        Color kD = Color(1.0f - metalness); //Metals do not have a diffuse part
        kD *= Color(1.0f) - F;//Only the transmitted light is diffused

        Color diffuse_part = kD * base_color / M_PI;
        Color specular_part = (F * D * G) / (4.0f * NoV * NoL);

        brdf_color = diffuse_part + specular_part;

        pdf = D * NoH / (4.0f * VoH);
    }

    m_out_stream << brdf_color << sycl::endl;
    m_out_stream << pdf << sycl::endl;
    return brdf_color;
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

bool RenderKernel::intersect_scene_bvh(Ray& ray, HitInfo& closest_hit_info) const
{
    closest_hit_info.t = -1;

    FlattenedBVH::Stack stack;
    stack.push(0);//Pushing the root of the BVH

    sycl::marray<float, BVHConstants::PLANES_COUNT> denoms;
    sycl::marray<float, BVHConstants::PLANES_COUNT> numers;

    for (int i = 0; i < BVHConstants::PLANES_COUNT; i++)
    {
        denoms[i] = dot(BVH_PLANE_NORMALS[i], ray.direction);
        numers[i] = dot(BVH_PLANE_NORMALS[i], Vector(ray.origin));
    }

    float closest_intersection_distance = -1;
    while (!stack.empty())
    {
        int node_index = stack.pop();
        const FlattenedBVH::FlattenedNode& node = m_bvh_nodes[node_index];

        if (node.intersect_volume(denoms, numers))
        {
            if (node.is_leaf)
            {
                for (int i = 0; i < node.nb_triangles; i++)
                {
                    int triangle_index = node.triangles_indices[i];

                    HitInfo local_hit_info;
                    if (m_triangle_buffer_access[triangle_index].intersect(ray, local_hit_info))
                    {
                        if (closest_intersection_distance > local_hit_info.t || closest_intersection_distance == -1)
                        {
                            closest_intersection_distance = local_hit_info.t;
                            closest_hit_info = local_hit_info;
                            closest_hit_info.triangle_index = triangle_index;
                        }
                    }
                }
            }
            else
                for (int i = 0; i < 8; i++)
                    stack.push(node.children[i]);
        }
    }

    return closest_hit_info.t > -1;
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
    intersect_scene_bvh(ray, hit_info);
    if (hit_info.t + 1.0e-4f < t_max)
        //There is something in between the light and the origin of the ray
        return true;
    else
        return false;
}
