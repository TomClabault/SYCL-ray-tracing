#include "render_kernel.h"

#include "triangle.h"

void branchlessONB(const Vector& n, Vector& b1, Vector& b2)
{
    float sign = std::copysign(1.0f, n.z);
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

    float phi = 2.0f * (float)M_PI * rand_1;
    float root = std::sqrt(1.0f - rand_2 * rand_2);

    pdf = 1.0f / (2.0f * (float)M_PI);

    //Generating a random direction in a local space with Z as the Up vector
    Vector random_dir_local_space(std::cos(phi) * root, std::sin(phi) * root, rand_2);
    return rotate_vector_around_normal(normal, random_dir_local_space);
}

Vector RenderKernel::cosine_weighted_direction_around_normal(const Vector& normal, float& pdf, xorshift32_generator& random_number_generator) const
{
    float rand_1 = random_number_generator();
    float rand_2 = random_number_generator();

    float sqrt_rand_2 = std::sqrt(rand_2);
    float phi = 2.0f * (float)M_PI * rand_1;
    float cos_theta = sqrt_rand_2;
    float sin_theta = std::sqrt(std::max(0.0f, 1.0f - cos_theta * cos_theta));

    pdf = sqrt_rand_2 / (float)M_PI;

    //Generating a random direction in a local space with Z as the Up vector
    Vector random_dir_local_space = Vector(std::cos(phi) * sin_theta, std::sin(phi) * sin_theta, sqrt_rand_2);
    return rotate_vector_around_normal(normal, random_dir_local_space);
}

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
    xorshift32_generator random_number_generator(31 + x * y * SAMPLES_PER_KERNEL);
    //Generating some numbers to make sure the generators of each thread spread apart
    //If not doing this, the generator shows clear artifacts until it has generated
    //a few numbers
    random_number_generator();
    random_number_generator();
    random_number_generator();

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
                bool intersection_found = INTERSECT_SCENE(ray, closest_hit_info);

                if (intersection_found)
                {
                    int material_index = m_materials_indices_buffer[closest_hit_info.primitive_index];
                    SimpleMaterial material = m_materials_buffer_access[material_index];

                    // --------------------------------------------------- //
                    // ---------- Direct lighting light sources ---------- //
                    // --------------------------------------------------- //
                    Color light_sources_radiance = sample_light_sources(ray, closest_hit_info, material, random_number_generator);
                    Color env_map_radiance = sample_environment_map(ray, closest_hit_info, material, random_number_generator);

                    // --------------------------------------- //
                    // ---------- Indirect lighting ---------- //
                    // --------------------------------------- //

                    float pdf;
                    Vector random_bounce_direction;// = cosine_weighted_direction_around_normal(closest_hit_info.normal_at_intersection, pdf, random_number_generator);
                    Color brdf = cook_torrance_brdf_importance_sample(material, -ray.direction, closest_hit_info.normal_at_intersection, random_bounce_direction, pdf, random_number_generator);
                    //Color brdf = cook_torrance_brdf(material, random_bounce_direction, -ray.direction, closest_hit_info.normal_at_intersection);
                    throughput *= brdf * std::max(0.0f, dot(random_bounce_direction, closest_hit_info.normal_at_intersection));
                    
                    if (bounce == 0)
                        sample_color += material.emission;
                    sample_color += (light_sources_radiance + env_map_radiance) * throughput;

                    //throughput /= pdf;

                    Point new_ray_origin = closest_hit_info.inter_point + closest_hit_info.normal_at_intersection * 1.0e-4f;
                    ray = Ray(new_ray_origin, random_bounce_direction);
                    next_ray_state = RayState::BOUNCE;
                }
                else
                    next_ray_state = RayState::MISSED;
            }
            else if (next_ray_state == MISSED)
            {
                /*
                Color skysphere_color = sample_environment_map_from_direction(ray.direction);

                sample_color += skysphere_color * throughput;
                */

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

    const float gamma = 2.2f;
    const float exposure = 0.75f;
    Color hdrColor = m_frame_buffer_access[y * m_width + x];

    //Exposure tone mapping
    Color tone_mapped = Color(1.0f, 1.0f, 1.0f) - exp(-hdrColor * exposure);
    // Gamma correction
    Color gamma_corrected = pow(tone_mapped, 1.0f / gamma);

    m_frame_buffer_access[y * m_width + x] = gamma_corrected;
}

#include <omp.h>
void RenderKernel::render()
{
    std::atomic<int> lines_completed = 0;

#pragma omp parallel for schedule(dynamic)
    for (int y = 0; y < m_frame_buffer_access.height(); y++)
    {
        for (int x = 0; x < m_frame_buffer_access.width(); x++)
            ray_trace_pixel(x, y);

        lines_completed++;

        if (omp_get_thread_num() == 0)
        {
            if (lines_completed % (m_frame_buffer_access.height() / 25))
                std::cout << lines_completed / (float)m_frame_buffer_access.height() * 100 << "%" << std::endl;
        }
    }
}

Color RenderKernel::lambertian_brdf(const SimpleMaterial& material, const Vector& to_light_direction, const Vector& view_direction, const Vector& surface_normal) const
{
    return material.diffuse / (float)M_PI;
}

Color fresnel_schlick(Color F0, float NoV)
{
    return F0 + (Color(1.0f) - F0) * std::pow((1.0f - NoV), 5.0f);
}

float GGX_normal_distribution(float alpha, float NoH)
{
    float alpha2 = alpha * alpha;
    float NoH2 = NoH * NoH;
    float b = (NoH2 * (alpha2 - 1.0f) + 1.0f);
    return alpha2 * (1.0f / (float)M_PI) / std::max(b * b, 0.0001f);
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

inline Color RenderKernel::cook_torrance_brdf(const SimpleMaterial& material, const Vector& to_light_direction, const Vector& view_direction, const Vector& surface_normal) const
{
    Color brdf_color = Color(0.0f, 0.0f, 0.0f);
    Color base_color = material.diffuse;

    Vector halfway_vector = normalize(view_direction + to_light_direction);

    float NoV = std::max(0.0f, dot(surface_normal, view_direction));
    float NoL = std::max(0.0f, dot(surface_normal, to_light_direction));
    float NoH = std::max(0.0f, dot(surface_normal, halfway_vector));
    float VoH = std::max(0.0f, dot(halfway_vector, view_direction));

    if (NoV > 0.0f && NoL > 0.0f && NoH > 0.0f)
    {
        float metalness = material.metalness;
        float roughness = material.roughness;

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

        Color diffuse_part = kD * base_color / (float)M_PI;
        Color specular_part = (F * D * G) / (4.0f * NoV * NoL);

        brdf_color = diffuse_part + specular_part;
    }

    return brdf_color;
}

Color fresnelSchlick(Color F0, float cosTheta) {
  return F0 + (Color(1.0f) - F0) * std::pow(1.0f - cosTheta, 5.0f);
}

float D_GGX(float alpha, float NoH) {
  float alpha2 = alpha * alpha;
  float NoH2 = NoH * NoH;
  float b = (NoH2 * (alpha2 - 1.0f) + 1.0f);
  return alpha2 * (1.0f / (float)M_PI) / (b * b);
}

float G1_GGX_Schlick(float NoV, float alpha) {
  float k = alpha / 2.0f;
  return std::max(NoV, 0.001f) / (NoV * (1.0f - k) + k);
}

float G_Smith(float alpha, float NoV, float NoL) {
  return G1_GGX_Schlick(NoL, alpha) * G1_GGX_Schlick(NoV, alpha);
}

inline Color RenderKernel::cook_torrance_brdf_importance_sample(const SimpleMaterial& material, const Vector& view_direction, const Vector& surface_normal, Vector& output_direction, float& pdf, xorshift32_generator& random_number_generator) const
{
    float metalness = material.metalness;
    float roughness = material.roughness;
    float alpha = roughness * roughness;

    float rand1 = random_number_generator();
    float rand2 = random_number_generator();

    float phi = 2.0f * (float)M_PI * rand1;
    float theta = std::acos((1.0f - rand2) / (rand2 * (alpha * alpha - 1.0f) + 1.0f));
    float sin_theta = std::sin(theta);

    Vector microfacet_normal_local_space = Vector(std::cos(phi) * sin_theta, std::sin(phi) * sin_theta, std::cos(theta));
    Vector microfacet_normal = rotate_vector_around_normal(surface_normal, microfacet_normal_local_space);
    if (dot(microfacet_normal, surface_normal) < 0.0f)
        //The microfacet normal that we sampled was under the surface, it can happen
        return Color(0.0f);
    Vector to_light_direction = 2.0f * dot(microfacet_normal, view_direction) * microfacet_normal - view_direction;
    Vector halfway_vector = microfacet_normal;
    output_direction = to_light_direction;

    Color brdf_color = Color(0.0f, 0.0f, 0.0f);
    Color base_color = material.diffuse;

    float NoV = std::max(0.0f, dot(surface_normal, view_direction));
    float NoL = std::max(0.0f, dot(surface_normal, to_light_direction));
    float NoH = std::max(0.0f, dot(surface_normal, halfway_vector));
    float VoH = std::max(0.0f, dot(halfway_vector, view_direction));

    if (NoV > 0.0f && NoL > 0.0f && NoH > 0.0f)
    {
        /////////// Cook Torrance BRDF //////////
        Color F;
        float D, G;

        //TODO check metalness parce que avec metalness = 1.0f, on a quand même des reflets blanc sur le mur vert
        //On devrait pas avoir des reflets verts metallic si on a metalness a 1.0 ?
        //F0 = 0.04 for dielectrics, 1.0 for metals (approximation)
        Color F0 = Color(0.04f * (1.0f - metalness)) + metalness * base_color;

        //GGX Distribution function
        F = fresnel_schlick(F0, VoH);
        D = GGX_normal_distribution(alpha, NoH);
        G = GGX_smith_masking_shadowing(alpha, NoV, NoL);

        Color kD = Color(1.0f - metalness); //Metals do not have a diffuse part
        kD *= Color(1.0f) - F;//Only the transmitted light is diffused

        Color diffuse_part = kD * base_color / (float)M_PI;
        Color specular_part = (F * D * G) / (4.0f * NoV * NoL);

        pdf = D * NoH / (4.0f * VoH);

        brdf_color = diffuse_part + specular_part / pdf;
    }

    return brdf_color;
}

bool RenderKernel::intersect_scene(const Ray& ray, HitInfo& closest_hit_info) const
{
    closest_hit_info.t = -1.0f;

    for (int i = 0; i < m_triangle_buffer_access.size(); i++)
    {
        const Triangle& triangle = m_triangle_buffer_access[i];

        HitInfo hit_info;
        if(triangle.intersect(ray, hit_info))
        {
            if (hit_info.t < closest_hit_info.t || closest_hit_info.t == -1.0f)
            {
                closest_hit_info = hit_info;
                closest_hit_info.primitive_index = i;
            }
        }
    }

    for (int i = 0; i < m_sphere_buffer.size(); i++)
    {
        const Sphere& sphere = m_sphere_buffer[i];

        HitInfo hit_info;
        if (sphere.intersect(ray, hit_info))
            if (hit_info.t < closest_hit_info.t || closest_hit_info.t == -1.0f)
                closest_hit_info = hit_info;
    }

    return closest_hit_info.t > 0.0f;
}

/*
* Flat BVH intersection for the GPU
*/
//inline bool RenderKernel::intersect_scene_bvh(const Ray& ray, HitInfo& closest_hit_info) const
//{
//    closest_hit_info.t = -1.0f;
//
//    FlattenedBVH::Stack stack;
//    stack.push(0);//Pushing the root of the BVH
//
//    std::array<float, BVHConstants::PLANES_COUNT> denoms;
//    std::array<float, BVHConstants::PLANES_COUNT> numers;
//
//    for (int i = 0; i < BVHConstants::PLANES_COUNT; i++)
//    {
//        denoms[i] = dot(BoundingVolume::PLANE_NORMALS[i], ray.direction);
//        numers[i] = dot(BoundingVolume::PLANE_NORMALS[i], Vector(ray.origin));
//    }
//
//    float closest_intersection_distance = -1;
//    while (!stack.empty())
//    {
//        int node_index = stack.pop();
//        const FlattenedBVH::FlattenedNode& node = m_bvh_nodes[node_index];
//
//        if (node.intersect_volume(denoms, numers))
//        {
//            if (node.is_leaf)
//            {
//                for (int i = 0; i < node.nb_triangles; i++)
//                {
//                    int triangle_index = node.triangles_indices[i];
//
//                    HitInfo local_hit_info;
//                    if (m_triangle_buffer_access[triangle_index].intersect(ray, local_hit_info))
//                    {
//                        if (closest_intersection_distance > local_hit_info.t || closest_intersection_distance == -1)
//                        {
//                            closest_intersection_distance = local_hit_info.t;
//                            closest_hit_info = local_hit_info;
//                            closest_hit_info.material_index = m_materials_indices_buffer[triangle_index];
//                        }
//                    }
//                }
//            }
//            else
//            {
//                stack.push(node.children[0]);
//                stack.push(node.children[1]);
//                stack.push(node.children[2]);
//                stack.push(node.children[3]);
//                stack.push(node.children[4]);
//                stack.push(node.children[5]);
//                stack.push(node.children[6]);
//                stack.push(node.children[7]);
//            }
//        }
//    }
//
//    for (int i = 0; i < m_sphere_buffer.size(); i++)
//    {
//        const Sphere& sphere = m_sphere_buffer[i];
//        HitInfo hit_info;
//        if (sphere.intersect(ray, hit_info))
//        {
//            if (hit_info.t < closest_intersection_distance || closest_intersection_distance == -1.0f)
//            {
//                closest_intersection_distance = hit_info.t;
//                closest_hit_info = hit_info;
//            }
//        }
//    }
//
//    return closest_hit_info.t > -1.0f;
//}

inline bool RenderKernel::intersect_scene_bvh(const Ray& ray, HitInfo& closest_hit_info) const
{
    closest_hit_info.t == -1.0f;

    m_bvh.intersect(ray, closest_hit_info);

    for (int i = 0; i < m_sphere_buffer.size(); i++)
    {
        const Sphere& sphere = m_sphere_buffer[i];

        HitInfo hit_info;
        if (sphere.intersect(ray, hit_info))
            if (hit_info.t < closest_hit_info.t || closest_hit_info.t == -1.0f)
                closest_hit_info = hit_info;
    }

    return closest_hit_info.t > 0.0f;
}

inline bool RenderKernel::INTERSECT_SCENE(const Ray& ray, HitInfo& hit_info) const
{
#if USE_BVH
    return intersect_scene_bvh(ray, hit_info);
#else
    return intersect_scene(ray, hit_info);
#endif
}

float power_heuristic(float pdf_a, float pdf_b)
{
    float pdf_a_squared = pdf_a * pdf_a;

    return pdf_a_squared / (pdf_a_squared + pdf_b * pdf_b);
}

Color RenderKernel::sample_environment_map_from_direction(const Vector& direction) const
{
    float u, v;
    u = 0.5f + std::atan2(direction.z, direction.x) / (2.0f * (float)M_PI);
    v = 0.5f + std::asin(direction.y) / (float)M_PI;

    int x = u * m_skysphere.width();
    int y = m_skysphere.height() - v * m_skysphere.height();

    return m_skysphere[y * m_skysphere.width() + x];
}

Color RenderKernel::sample_environment_map(const Ray& ray, const HitInfo& closest_hit_info, const SimpleMaterial& material, xorshift32_generator& random_number_generator) const
{
    Vector sampled_direction;

    int env_map_bin_area;
    int random_pixel_coord_x, random_pixel_coord_y;
    float sin_theta;

    do
    {
        int random_bin_index = m_env_map_bins.size() * random_number_generator();
        ImageBin env_map_bin = m_env_map_bins[random_bin_index];

        int bin_width = env_map_bin.x1 - env_map_bin.x0;
        int bin_height = env_map_bin.y1 - env_map_bin.y0;
        env_map_bin_area = bin_width * bin_height;

        random_pixel_coord_x = random_number_generator() * bin_width + env_map_bin.x0;
        random_pixel_coord_y = random_number_generator() * bin_height + env_map_bin.y0;

        //Converting from pixel coordinates to world direction
        float u = (float)random_pixel_coord_x / m_skysphere.width();
        float v = (float)random_pixel_coord_y / m_skysphere.height();

        float phi = u * 2.0f * M_PI;
        float theta = v * M_PI;

        sin_theta = std::sin(theta);
        float cos_theta = std::sqrt(1 - sin_theta * sin_theta);

        // Convert to cartesian coordinates
        sampled_direction = Vector(std::cos(phi) * sin_theta, std::sin(phi) * sin_theta, cos_theta);

        //While we sampled a direction behind the object (our objects are opaque)
        //so we're not going to see the environment map if it's below the surface
        //of the object
    } while (dot(sampled_direction, closest_hit_info.normal_at_intersection) < 0);

    if (evaluate_shadow_ray(Ray(ray.origin + closest_hit_info.normal_at_intersection * 1.0e-5f, sampled_direction), (float)-INFINITY))
        return Color(0.0f);
    else
    {
        if (sin_theta == 0.0f)
            return Color(0.0f);
        else
        {
            float pixel_pdf = (float)(m_skysphere.width() * m_skysphere.height()) / (m_env_map_bins.size() * env_map_bin_area);
            //Conversion to solid angle
            float skysphere_pdf = pixel_pdf / (2.0f * M_PI * M_PI * sin_theta);

            Vector brdf_sample_direction;
            float brdf_pdf;
            Color brdf = cook_torrance_brdf_importance_sample(material, -ray.direction, closest_hit_info.normal_at_intersection, brdf_sample_direction, brdf_pdf, random_number_generator);

            //Incoming radiance in the direction importance sampled from the env map pdf
            Color skysphere_sample_radiance = m_skysphere(random_pixel_coord_x, random_pixel_coord_y);

            //Incoming radiance in the direction importance sampled from the brdf pdf
            Color brdf_sample_radiance = sample_environment_map_from_direction(brdf_sample_direction);

            //Multiple Importance Sampling between the two samples for the environment map lighting
            float mis_weight = power_heuristic(brdf_pdf, skysphere_pdf);
            //return brdf_sample_radiance * mis_weight * (skysphere_sample_radiance * brdf) / skysphere_pdf;
            return (brdf_sample_radiance + (skysphere_sample_radiance * brdf) / skysphere_pdf) / 2;
        }
    }
}

Color RenderKernel::sample_light_sources(const Ray& ray, const HitInfo& closest_hit_info, const SimpleMaterial& material, xorshift32_generator& random_number_generator) const
{
    Color radiance = Color(0.0f, 0.0f, 0.0f);

    if (m_emissive_triangle_indices_buffer.size() > 0)
    {
        float light_sample_pdf;
        LightSourceInformation light_source_info;
        Point random_light_point = sample_random_point_on_lights(random_number_generator, light_sample_pdf, light_source_info);
        Point shadow_ray_origin = closest_hit_info.inter_point + closest_hit_info.normal_at_intersection * 1.0e-4f;
        Vector shadow_ray_direction = random_light_point - shadow_ray_origin;
        float distance_to_light = length(shadow_ray_direction);
        Vector shadow_ray_direction_normalized = normalize(shadow_ray_direction);

        Ray shadow_ray(shadow_ray_origin, shadow_ray_direction_normalized);

        bool in_shadow = evaluate_shadow_ray(shadow_ray, distance_to_light);

        if (!in_shadow)
        {
            const SimpleMaterial& emissive_triangle_material = m_materials_buffer_access[m_materials_indices_buffer[light_source_info.emissive_triangle_index]];

            radiance = emissive_triangle_material.emission;
            //Cosine angle on the illuminated surface
            radiance *= std::max(dot(closest_hit_info.normal_at_intersection, shadow_ray_direction_normalized), 0.0f);
            //Cosine angle on the light surface
            radiance *= std::max(dot(light_source_info.light_source_normal, -shadow_ray_direction_normalized), 0.0f);
            //Falloff of the light intensity with the distance squared
            radiance /= distance_to_light * distance_to_light;
            //PDF: Probability of having chosen this point on this exact light source
            radiance /= light_sample_pdf;
            //BRDF of the illuminated surface
            radiance *= cook_torrance_brdf(material, shadow_ray.direction, -ray.direction, closest_hit_info.normal_at_intersection);
        }
    }

    return radiance;
}

inline Point RenderKernel::sample_random_point_on_lights(xorshift32_generator& random_number_generator, float& pdf, LightSourceInformation& light_info) const
{
    light_info.emissive_triangle_index = random_number_generator() * m_emissive_triangle_indices_buffer.size();
    light_info.emissive_triangle_index = m_emissive_triangle_indices_buffer[light_info.emissive_triangle_index];
    Triangle random_emissive_triangle = m_triangle_buffer_access[light_info.emissive_triangle_index];

    float rand_1 = random_number_generator();
    float rand_2 = random_number_generator();

    float sqrt_r1 = std::sqrt(rand_1);
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

bool RenderKernel::evaluate_shadow_ray(const Ray& ray, float t_max) const
{
    HitInfo hit_info;
    bool inter_found = INTERSECT_SCENE(ray, hit_info);

    if (inter_found)
    {
        if (hit_info.t + 1.0e-4f < t_max)
            //There is something in between the light and the origin of the ray
            return true;
        else
            return false;
    }
}
