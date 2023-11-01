#ifndef RENDER_KERNEL_H
#define RENDER_KERNEL_H

#include "camera.h"
#include "color.h"
#include "triangle.h"

#include <sycl/sycl.hpp>

#define SAMPLES_PER_KERNEL 64
#define MAX_BOUNCES 2

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
                 sycl::accessor<Color, 1, sycl::access::mode::write, sycl::access::target::device> frame_buffer_accessor,
                 sycl::accessor<Triangle, 1, sycl::access::mode::read, sycl::access::target::device> triangle_buffer_accessor) :
        m_width(width),
        m_height(height), 
        m_frame_buffer_access(frame_buffer_accessor),
        m_triangle_buffer_access(triangle_buffer_accessor) {}

    void set_camera(Camera camera) { m_camera = camera; }

    void operator()(const sycl::nd_item<2>& coordinates) const
    {
        int x = coordinates.get_global_id(0);
        int y = coordinates.get_global_id(1);

        ray_trace_pixel(x, y);
    }

    SYCL_EXTERNAL Ray get_camera_ray(float x, float y) const;

    SYCL_EXTERNAL void ray_trace_pixel(int x, int y) const;

    bool intersect_scene(const Ray ray, HitInfo* closest_hit_info) const;
    bool evaluate_shadow_ray(Ray ray, float t_max) const;

private:
    int m_width, m_height;

    sycl::accessor<Color, 1, sycl::access::mode::write, sycl::access::target::device> m_frame_buffer_access;

    sycl::accessor<Triangle, 1, sycl::access::mode::read, sycl::access::target::device> m_triangle_buffer_access;

    Camera m_camera;
};

#endif

