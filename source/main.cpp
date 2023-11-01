#include <iostream>
#include <sycl/sycl.hpp>
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include "camera.h"
#include "render_kernel.h"
#include "triangle.h"

inline float clamp(const float x, const float min, const float max)
{
    if (x < min) return min;
    else if (x > max) return max;
    else return x;
}

bool write_image_png(std::vector<Color>& image, int width, int height, const char* filename)
{
    if (image.size() == 0)
        return false;

    std::vector<unsigned char> tmp(width * height * 4);
    for (unsigned i = 0, offset = 0; i < image.size(); i++, offset += 4)
    {
        Color pixel = image[i] * 255;
        tmp[offset] = clamp(pixel.r, 0, 255);
        tmp[offset + 1] = clamp(pixel.g, 0, 255);
        tmp[offset + 2] = clamp(pixel.b, 0, 255);
        tmp[offset + 3] = clamp(pixel.a, 0, 255);
    }

    stbi_flip_vertically_on_write(true);
    return stbi_write_png(filename, width, height, 4, tmp.data(), width * 4) != 0;
}

int main(int argc, char* argv[])
{
    const int width = 1280;
    const int height = 720;

    sycl::queue queue {sycl::cpu_selector_v};
    std::cout << "Using " << queue.get_device().get_info<sycl::info::device::name>() << std::endl;

    std::vector<Color> image(width * height);

    sycl::buffer<Color> image_buffer(image.data(), width * height);
    std::vector<Triangle> triangles{ Triangle(Point(-1, 0, -2), Point(1, 0, -2), Point(0.5, 0.5, -2)) };
    sycl::buffer<Triangle> triangle_buffer(triangles.data(), triangles.size());

    std::cout << "[" << width << "x" << height << "]: " << SAMPLES_PER_KERNEL << " samples" << std::endl << std::endl;

    queue.submit([&](sycl::handler& handler) {
        auto image_buffer_access = image_buffer.get_access<sycl::access::mode::write>(handler);
        auto triangle_buffer_access = triangle_buffer.get_access<sycl::access::mode::read>(handler);

        const auto global_range = sycl::range<2>(width, height);
        const auto local_range = sycl::range<2>(TILE_SIZE_X, TILE_SIZE_Y);
        const auto coordinates_indices = sycl::nd_range<2>(global_range, local_range);

        auto render_kernel = RenderKernel(width, height,
            image_buffer_access,
            triangle_buffer_access);
        render_kernel.set_camera(Camera(45, Translation(0, 1, 3.5)));

        handler.parallel_for(coordinates_indices, render_kernel);
    }).wait();

    image_buffer.get_host_access();

    write_image_png(image, width, height, "../TP_RT_output.png");

    return 0;
}
