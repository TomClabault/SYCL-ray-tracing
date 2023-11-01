#include <iostream>
#include <sycl/sycl.hpp>
#include <chrono>
#include <cmath>
#include <vector>

#include <rapidobj.hpp>

#include <stb_image_write.h>

#include "camera.h"
#include "render_kernel.h"
#include "triangle.h"

void write_image_png(std::vector<Color> image, const std::string& filepath)
{

}

int main(int argc, char* argv[])
{
    const int width = 1280;
    const int height = 720;
    sycl::queue queue {sycl::cpu_selector_v};
    std::cout << "Using " << queue.get_device().get_info<sycl::info::device::name>() << std::endl;

    std::vector<Color> image(width * height);

    sycl::buffer<Color> image_buffer(image.data(), width * height);
    std::vector<Triangle> triangles;
    triangles.push_back(Triangle(Point(-1, 0, -2), Point(1, 0, -2), Point(0.5, 0.5, -2)));
    sycl::buffer<Triangle> triangle_buffer(triangles.data(), triangles.size());

    std::cout << "[" << width << "x" << height << "]: " << RENDER_KERNEL_ITERATIONS * SAMPLES_PER_KERNEL << " samples" << std::endl << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    for (int kernel_iteration = 0; kernel_iteration < RENDER_KERNEL_ITERATIONS; kernel_iteration++)
    {
        try
        {
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
        }
        catch (sycl::exception e)
        {
            std::cout << "Kernel exception: ";
            std::cout << e.what() << std::endl;
            std::exit(-1);
        }

        std::cout << (float)(kernel_iteration + 1) / RENDER_KERNEL_ITERATIONS * 100.0f << "%" << std::endl;
    }

    queue.wait();

    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << "ms" << std::endl;

    image_buffer.get_host_access();

    write_image_png(image, "../TP_RT_output.png");

    return 0;
}
