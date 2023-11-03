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

#define SAMPLES_PER_KERNEL 2048
#define MAX_BOUNCES 5

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
                 Image& skysphere) : 
        m_width(width), m_height(height),
        m_frame_buffer_access(image_buffer),
        m_triangle_buffer_access(triangle_buffer_accessor),
        m_materials_buffer_access(materials_buffer_accessor),
        m_emissive_triangle_indices_buffer(emissive_triangle_indices_buffer_accessor),
        m_materials_indices_buffer(materials_indices_buffer_accessor),
        m_sphere_buffer(analytic_spheres_buffer),
        m_bvh(bvh),
        m_skysphere(skysphere) {}

    void set_camera(Camera camera) { m_camera = camera; }

    Ray get_camera_ray(float x, float y) const;

    Vector rotate_vector_around_normal(const Vector& normal, const Vector& random_dir_local_space) const;
    Vector uniform_direction_around_normal(const Vector& normal, float& pdf, xorshift32_generator& random_number_generator) const;
    Vector cosine_weighted_direction_around_normal(const Vector& normal, float& pdf, xorshift32_generator& random_number_generator) const;

    void ray_trace_pixel(int x, int y) const;
    void render();

    Color lambertian_brdf(const SimpleMaterial& material, const Vector& to_light_direction, const Vector& view_direction, const Vector& surface_normal) const;
    Color cook_torrance_brdf(const SimpleMaterial& material, const Vector& to_light_direction, const Vector& view_direction, const Vector& surface_normal) const;
    Color cook_torrance_brdf_importance_sample(const SimpleMaterial& material, const Vector& view_direction, const Vector& surface_normal, Vector& output_direction, xorshift32_generator& random_number_generator) const;

    bool intersect_scene(const Ray& ray, HitInfo& closest_hit_info) const;
    bool intersect_scene_bvh(const Ray& ray, HitInfo& closest_hit_info) const;
    bool INTERSECT_SCENE(const Ray& ray, HitInfo& hit_info)const ;
    Point sample_random_point_on_lights(xorshift32_generator& random_number_generator, float& pdf, LightSourceInformation& light_info) const;
    bool evaluate_shadow_ray(const Ray& ray, float t_max) const;

private:
    int m_width, m_height;
    int m_kernel_iteration;

    Image& m_frame_buffer_access;

    std::vector<Triangle> m_triangle_buffer_access;
    std::vector<SimpleMaterial> m_materials_buffer_access;
    std::vector<int> m_emissive_triangle_indices_buffer;
    std::vector<int> m_materials_indices_buffer;

    std::vector<Sphere> m_sphere_buffer;

    BVH& m_bvh;

    Image& m_skysphere;

    Camera m_camera;
};

#endif

