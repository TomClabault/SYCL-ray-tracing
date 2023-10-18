#include <iostream>
#include <sycl/sycl.hpp>
#include <chrono>
#include <cmath>

#include <rapidobj/rapidobj.hpp>

#include "camera.h"
#include "image_io.h"
#include "render_kernel.h"
#include "simple_material.h"
#include "tests.h"
#include "triangle.h"

#include "xorshift.h"

int dichotomie(std::vector<float> bins, float random)
{
    return 0;
}

int main(int argc, char* argv[])
{
    regression_tests();

    std::vector<float> bins {4, 5, 8, 2, 1};
    std::vector<float> bins2 (bins.size());
    std::vector<int> results(bins.size());
    float sum = 0;
    std::for_each(bins.begin(), bins.end(), [&sum] (float element) {
        sum += element;
    });

    float cumul = 0.0f;
    for (int i = 0; i < bins.size(); i++)
    {
        bins2[i] = bins[i] / sum + cumul;
        cumul += bins[i] / sum;
    }

    xorshift32_generator generator(58);


    for (int j = 0; j<  1000000; j++)
    {
        float random = generator();

        results[dichotomie(bins2, random)]++;

        /*
        for (int i = 0; i < bins.size(); i++)
        {
            if (random < bins2[i])
            {
                results[i]++;

                break;
            }
        }
        */
    }

    //return 0;

    const int width = 1280;
    const int height = 720;

    sycl::queue queue {sycl::gpu_selector_v};
    std::cout << "Using " << queue.get_device().get_info<sycl::info::device::name>() << std::endl;

    Image image(width, height);

    rapidobj::Result parsed_obj = rapidobj::ParseFile("../data/cornell.obj", rapidobj::MaterialLibrary::Default());
    //rapidobj::Result parsed_obj = rapidobj::ParseFile("../data/test_triangle_area_sampling.obj", rapidobj::MaterialLibrary::Default());
    if (parsed_obj.error)
    {
        std::cout << "There was an error loading the OBJ file: " << parsed_obj.error.code.message() << std::endl;

        return -1;
    }
    rapidobj::Triangulate(parsed_obj);

    const rapidobj::Array<float>& positions = parsed_obj.attributes.positions;
    std::vector<Triangle> triangle_host_buffer;
    std::vector<int> emissive_triangle_indices_host_buffer;
    std::vector<int> materials_indices_host_buffer;
    for (rapidobj::Shape& shape : parsed_obj.shapes)
    {
        rapidobj::Mesh& mesh = shape.mesh;
        for (int i = 0; i < mesh.indices.size(); i += 3)
        {
            int index_0 = mesh.indices[i + 0].position_index;
            int index_1 = mesh.indices[i + 1].position_index;
            int index_2 = mesh.indices[i + 2].position_index;

            Point A = Point(positions[index_0 * 3 + 0], positions[index_0 * 3 + 1], positions[index_0 * 3 + 2]);
            Point B = Point(positions[index_1 * 3 + 0], positions[index_1 * 3 + 1], positions[index_1 * 3 + 2]);
            Point C = Point(positions[index_2 * 3 + 0], positions[index_2 * 3 + 1], positions[index_2 * 3 + 2]);

            Triangle triangle(A, B, C);
            triangle_host_buffer.push_back(triangle);

            int mesh_triangle_index = i / 3;
            int triangle_material_index = mesh.material_ids[mesh_triangle_index];
            materials_indices_host_buffer.push_back(triangle_material_index);

            rapidobj::Float3 emission = parsed_obj.materials[triangle_material_index].emission;
            if (emission[0] > 0 || emission[1] > 0 || emission[2] > 0)
            {
                //This is an emissive triangle
                //Using the buffer of all the triangles to get the global id of the emissive
                //triangle
                emissive_triangle_indices_host_buffer.push_back(triangle_host_buffer.size() - 1);
            }
        }
    }

    //Computing SimpleMaterials
    std::vector<SimpleMaterial> materials_host_buffer;
    for (const rapidobj::Material& material : parsed_obj.materials)
        materials_host_buffer.push_back(SimpleMaterial {Color(material.emission), Color(material.diffuse)});

    sycl::buffer<Color> image_buffer(image.color_data(), image.width() * image.height());
    sycl::buffer<Triangle> triangle_buffer(triangle_host_buffer.data(), triangle_host_buffer.size());
    sycl::buffer<SimpleMaterial> materials_buffer(materials_host_buffer.data(), materials_host_buffer.size());
    sycl::buffer<int> emissive_triangle_indices_buffer(emissive_triangle_indices_host_buffer.data(), emissive_triangle_indices_host_buffer.size());
    sycl::buffer<int> materials_indices_buffer(materials_indices_host_buffer.data(), materials_indices_host_buffer.size());

    std::cout << "[" << width << "x" << height << "]: " << RENDER_KERNEL_ITERATIONS * SAMPLES_PER_KERNEL << " samples" << std::endl << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < RENDER_KERNEL_ITERATIONS; i++)
    {
        queue.submit([&] (sycl::handler& handler) {
            auto image_buffer_access = image_buffer.get_access<sycl::access::mode::write, sycl::access::target::device>(handler);
            auto triangle_buffer_access = triangle_buffer.get_access<sycl::access::mode::read>(handler);
            auto materials_buffer_access = materials_buffer.get_access<sycl::access::mode::read>(handler);
            auto emissive_triangle_indices_buffer_access = emissive_triangle_indices_buffer.get_access<sycl::access::mode::read>(handler);
            auto materials_indices_buffer_access = materials_indices_buffer.get_access<sycl::access::mode::read>(handler);

            const auto global_range = sycl::range<2>(width, height);
            const auto local_range = sycl::range<2>(TILE_SIZE_X, TILE_SIZE_Y);
            const auto coordinates_indices = sycl::nd_range<2>(global_range, local_range);

            sycl::stream debug_out_stream(1024, 128, handler);

            auto render_kernel = RenderKernel(width, height, i,
                                              image_buffer_access,
                                              triangle_buffer_access,
                                              materials_buffer_access,
                                              emissive_triangle_indices_buffer_access,
                                              materials_indices_buffer_access,
                                              debug_out_stream);
            render_kernel.set_camera(Camera(45, Translation(0, 1, 3.5)));

            handler.parallel_for(coordinates_indices, render_kernel);
        }).wait();

        std::cout << (float)(i + 1) / RENDER_KERNEL_ITERATIONS * 100.0f << "%" << std::endl;
    }

    queue.wait();

    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << "ms" << std::endl;

    image_buffer.get_access<sycl::access::mode::read>();


    write_image_png(image, "../TP_RT_output.png");

    return 0;
}
