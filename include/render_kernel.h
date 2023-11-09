#ifndef RENDER_KERNEL_H
#define RENDER_KERNEL_H

#include "bvh.h"
#include "camera.h"
#include "color.h"
#include "image.h"
#include "simple_material.h"
#include "sphere.h"
#include "triangle.h"
#include "xorshift.h"

#define SAMPLES_PER_KERNEL 32
#define MAX_BOUNCES 15
#define USE_BVH 1

#define TILE_SIZE_X 8
#define TILE_SIZE_Y TILE_SIZE_X

struct LightSourceInformation
{
    int emissive_triangle_index = -1;
    Vector light_source_normal;
};

class RenderKernel
{
public:
    RenderKernel(int width, int height,
                 Image& image_buffer,
                 const std::vector<Triangle>& triangle_buffer_accessor,
                 const std::vector<SimpleMaterial>& materials_buffer_accessor,
                 const std::vector<int>& emissive_triangle_indices_buffer_accessor,
                 const std::vector<int>& materials_indices_buffer_accessor,
                 const std::vector<Sphere>& analytic_spheres_buffer,
                 BVH& bvh,
                 const Image& skysphere,
                 const std::vector<ImageBin>& env_map_bins) : 
        m_width(width), m_height(height),
        m_frame_buffer(image_buffer),
        m_triangle_buffer_access(triangle_buffer_accessor),
        m_materials_buffer_access(materials_buffer_accessor),
        m_emissive_triangle_indices_buffer(emissive_triangle_indices_buffer_accessor),
        m_materials_indices_buffer(materials_indices_buffer_accessor),
        m_sphere_buffer(analytic_spheres_buffer),
        m_bvh(bvh),
        m_skysphere(skysphere),
        m_env_map_bins(env_map_bins) {}

    void set_camera(Camera camera) { m_camera = camera; }

    Ray get_camera_ray(float x, float y) const;

    Vector rotate_vector_around_normal(const Vector& normal, const Vector& random_dir_local_space) const;
    Vector uniform_direction_around_normal(const Vector& normal, float& pdf, xorshift32_generator& random_number_generator) const;
    Vector cosine_weighted_direction_around_normal(const Vector& normal, float& pdf, xorshift32_generator& random_number_generator) const;

    void ray_trace_pixel(int x, int y) const;
    void render();

    Color lambertian_brdf(const SimpleMaterial& material, const Vector& to_light_direction, const Vector& view_direction, const Vector& surface_normal) const;
    float cook_torrance_brdf_pdf(const SimpleMaterial& material, const Vector& view_direction, const Vector& to_light_direction, const Vector& surface_normal) const;
    Color cook_torrance_brdf(const SimpleMaterial& material, const Vector& to_light_direction, const Vector& view_direction, const Vector& surface_normal) const;
    Color cook_torrance_brdf_importance_sample(const SimpleMaterial& material, const Vector& view_direction, const Vector& surface_normal, Vector& output_direction, float& pdf, xorshift32_generator& random_number_generator) const;

    bool intersect_scene(const Ray& ray, HitInfo& closest_hit_info) const;
    bool intersect_scene_bvh(const Ray& ray, HitInfo& closest_hit_info) const;
    bool INTERSECT_SCENE(const Ray& ray, HitInfo& hit_info)const ;

    Color sample_environment_map_from_direction(const Vector& direction) const;
    Color sample_environment_map(const Ray& ray, const HitInfo& closest_hit_info, const SimpleMaterial& material, xorshift32_generator& random_number_generator) const;
    Color sample_light_sources(const Ray& ray, const HitInfo& closest_hit_info, const SimpleMaterial& material, xorshift32_generator& random_number_generator) const;
    Point sample_random_point_on_lights(xorshift32_generator& random_number_generator, float& pdf, LightSourceInformation& light_info) const;
    bool evaluate_shadow_ray(const Ray& ray, float t_max) const;

private:
    int m_width, m_height;

    Image& m_frame_buffer;

    const std::vector<Triangle>& m_triangle_buffer_access;
    const std::vector<SimpleMaterial>& m_materials_buffer_access;
    const std::vector<int>& m_emissive_triangle_indices_buffer;
    const std::vector<int>& m_materials_indices_buffer;

    const std::vector<Sphere>& m_sphere_buffer;

    const BVH& m_bvh;

    const Image& m_skysphere;
    const std::vector<ImageBin>& m_env_map_bins;

    Camera m_camera;
};

#endif

