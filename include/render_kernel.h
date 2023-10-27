#ifndef RENDER_KERNEL_H
#define RENDER_KERNEL_H

#include "bvh.h"
#include "camera.h"
#include "color.h"
#include "simple_material.h"
#include "triangle.h"
#include "xorshift.h"

#include <sycl/sycl.hpp>

#define RENDER_KERNEL_ITERATIONS 8
#define SAMPLES_PER_KERNEL 64
#define MAX_BOUNCES 5

#define TILE_SIZE_X 8
#define TILE_SIZE_Y TILE_SIZE_X

struct LightSourceInformation
{
    int emissive_triangle_index;
    Vector light_source_normal;
};

class RenderKernel
{
public:
    RenderKernel(int width, int height, int kernel_iteration,
                 sycl::accessor<Color, 1, sycl::access::mode::write, sycl::access::target::device> frame_buffer_accessor,
                 sycl::accessor<Triangle, 1, sycl::access::mode::read, sycl::access::target::device> triangle_buffer_accessor,
                 sycl::accessor<SimpleMaterial, 1, sycl::access::mode::read, sycl::access::target::device> materials_buffer_accessor,
                 sycl::accessor<int, 1, sycl::access::mode::read, sycl::access::target::device> emissive_triangle_indices_buffer_accessor,
                 sycl::accessor<int, 1, sycl::access::mode::read, sycl::access::target::device> materials_indices_buffer_accessor,
                 sycl::accessor<FlattenedBVH::FlattenedNode, 1, sycl::access::mode::read, sycl::access::target::device> bvh_nodes,
                 sycl::accessor<sycl::float4, 2, sycl::access::mode::read, sycl::access::target::image> skysphere,
                 sycl::sampler skysphere_sampler,
                 sycl::stream debug_out_stream) :
        m_width(width), m_height(height), m_kernel_iteration(kernel_iteration),
        m_frame_buffer_access(frame_buffer_accessor),
        m_triangle_buffer_access(triangle_buffer_accessor),
        m_materials_buffer_access(materials_buffer_accessor),
        m_emissive_triangle_indices_buffer(emissive_triangle_indices_buffer_accessor),
        m_materials_indices_buffer(materials_indices_buffer_accessor),
        m_bvh_nodes(bvh_nodes),
        m_skysphere(skysphere),
        m_skysphere_width(skysphere.get_range()[1]),
        m_skysphere_height(skysphere.get_range()[0]),
        m_skysphere_sampler(skysphere_sampler),
        m_out_stream(debug_out_stream) {}

    void set_camera(Camera camera) { m_camera = camera; }

    void set_bvh_plane_normals(const sycl::accessor<Vector, 1, sycl::access::mode::read, sycl::access::target::constant_buffer>& plane_normals_accessor)
    {
        BVH_PLANE_NORMALS = plane_normals_accessor;
    }

    SYCL_EXTERNAL void operator()(const sycl::nd_item<2>& coordinates) const;
    SYCL_EXTERNAL Ray get_camera_ray(float x, float y) const;

    Vector rotate_vector_around_normal(const Vector& normal, const Vector& random_dir_local_space) const;
    Vector cosine_weighted_direction_around_normal(const Vector& normal, float& pdf, xorshift32_generator& random_number_generator) const;
    Vector uniform_direction_around_normal(const Vector& normal, float& pdf, xorshift32_generator& random_number_generator) const;

    SYCL_EXTERNAL void ray_trace_pixel(int x, int y) const;

    SYCL_EXTERNAL Color lambertian_brdf(const SimpleMaterial& material, const Vector& to_light_direction, const Vector& view_direction, const Vector& surface_normal) const;
    SYCL_EXTERNAL Color microfacet_brdf(const SimpleMaterial& material, const Vector& to_light_direction, const Vector& view_direction, const Vector& surface_normal) const;

    SYCL_EXTERNAL bool intersect_scene(Ray& ray, HitInfo& closest_hit_info) const;
    SYCL_EXTERNAL bool intersect_scene_bvh(Ray& ray, HitInfo& closest_hit_info) const;
    Point sample_random_point_on_lights(xorshift32_generator& random_number_generator, float& pdf, LightSourceInformation& light_info) const;
    SYCL_EXTERNAL bool evaluate_shadow_ray(Ray& ray, float t_max) const;

private:
    int m_width, m_height;
    int m_kernel_iteration;

    sycl::accessor<Color, 1, sycl::access::mode::write, sycl::access::target::device> m_frame_buffer_access;

    sycl::accessor<Triangle, 1, sycl::access::mode::read, sycl::access::target::device> m_triangle_buffer_access;
    sycl::accessor<SimpleMaterial, 1, sycl::access::mode::read, sycl::access::target::device> m_materials_buffer_access;
    sycl::accessor<int, 1, sycl::access::mode::read, sycl::access::target::device> m_emissive_triangle_indices_buffer;
    sycl::accessor<int, 1, sycl::access::mode::read, sycl::access::target::device> m_materials_indices_buffer;

    sycl::accessor<FlattenedBVH::FlattenedNode, 1, sycl::access::mode::read, sycl::access::target::device> m_bvh_nodes;
    sycl::accessor<Vector, 1, sycl::access::mode::read, sycl::access::target::constant_buffer> BVH_PLANE_NORMALS;

    int m_skysphere_width, m_skysphere_height;
    sycl::accessor<sycl::float4, 2, sycl::access::mode::read, sycl::access::target::image> m_skysphere;
    sycl::sampler m_skysphere_sampler;

    sycl::stream m_out_stream;

    Camera m_camera;
};

#endif

