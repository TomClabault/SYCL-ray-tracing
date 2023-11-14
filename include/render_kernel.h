#ifndef RENDER_KERNEL_H
#define RENDER_KERNEL_H

#include "camera.h"
#include "color.h"
#include "triangle.h"

#define SAMPLES_PER_KERNEL 64
#define MAX_BOUNCES 4

#define TILE_SIZE_X 8
#define TILE_SIZE_Y TILE_SIZE_X

#include <vector>

struct LightSourceInformation
{
    int emissive_triangle_index = -1;
    Vector light_source_normal;
};

class RenderKernel
{
public:
    RenderKernel(int width, int height,
                 std::vector<Color>& frame_buffer_accessor,
                 std::vector<Triangle>& triangle_buffer_accessor) :
        m_width(width),
        m_height(height), 
        m_frame_buffer_access(frame_buffer_accessor),
        m_triangle_buffer_access(triangle_buffer_accessor) {}

    void set_camera(Camera camera) { m_camera = camera; }

    void render()
    {
        for (int y = 0; y < m_height; y++)
        {
            for (int x = 0; x < m_width; x++)
                ray_trace_pixel(x, y);

            std::cout << (float)y / m_height * 100 << "%" << std::endl;
        }
    }

    Ray get_camera_ray(float x, float y) const;

    void ray_trace_pixel(int x, int y);

    bool intersect_scene(const Ray ray, HitInfo& closest_hit_info) const;
    bool evaluate_shadow_ray(Ray ray, float t_max) const;

private:
    int m_width, m_height;

    std::vector<Color>& m_frame_buffer_access;
    std::vector<Triangle>& m_triangle_buffer_access;

    Camera m_camera;
};

#endif

