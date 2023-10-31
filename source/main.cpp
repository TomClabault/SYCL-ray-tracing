#include <iostream>
#include <sycl/sycl.hpp>
#include <chrono>
#include <cmath>

#include <rapidobj.hpp>

#include <stb_image_write.h>

#include "bvh.h"
#include "camera.h"
#include "image_io.h"
#include "render_kernel.h"
#include "simple_material.h"
#include "sphere.h"
#include "tests.h"
#include "triangle.h"
#include "utils.h"

#include "xorshift.h"

Sphere add_sphere_to_scene(ParsedOBJ& parsed_obj, const Point& center, float radius, const SimpleMaterial& material)
{
    int material_index = parsed_obj.materials.size();
    parsed_obj.materials.push_back(material);

    parsed_obj.material_indices.push_back(material_index);

    Sphere sphere(center, radius, material_index);

    return sphere;
}

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

    const int width = 1280;
    const int height = 720;
    sycl::queue queue {sycl::cpu_selector_v};
    std::cout << "Using " << queue.get_device().get_info<sycl::info::device::name>() << std::endl;

    Image image(width, height);

    ParsedOBJ parsed_obj = Utils::parse_obj("../SYCL-ray-tracing/data/OBJs/cornell_pbr.obj");

    /*Sphere sphere = add_sphere_to_scene(parsed_obj, Point(0.3275, 0.7, 0.3725), 0.2, SimpleMaterial {Color(0.0f), Color(1.0f, 0.71, 0.29), 1.0f, 0.4f});
    std::vector<Sphere> spheres = { sphere };*/

    BVH bvh(&parsed_obj.triangles);
    FlattenedBVH flat_bvh = bvh.flatten();

    sycl::buffer<Color> image_buffer(image.color_data(), image.width() * image.height());
    sycl::buffer<Triangle> triangle_buffer(parsed_obj.triangles.data(), parsed_obj.triangles.size());
    sycl::buffer<SimpleMaterial> materials_buffer(parsed_obj.materials.data(), parsed_obj.materials.size());
    sycl::buffer<int> emissive_triangle_indices_buffer(parsed_obj.emissive_triangle_indices.data(), parsed_obj.emissive_triangle_indices.size());
    sycl::buffer<int> materials_indices_buffer(parsed_obj.material_indices.data(), parsed_obj.material_indices.size());
    //sycl::buffer<Sphere> sphere_buffer(spheres.data(), spheres.size());
    sycl::buffer<FlattenedBVH::FlattenedNode> bvh_nodes_buffer(flat_bvh.get_nodes().data(), flat_bvh.get_nodes().size());
    sycl::buffer<Vector> bvh_plane_normals_buffer(BoundingVolume::PLANE_NORMALS, BVHConstants::PLANES_COUNT);

    //int skysphere_width, skysphere_height;
    //std::vector<sycl::float4> skysphere_data = Utils::read_image_float("../SYCL-ray-tracing/data/Skyspheres/evening_road_01_puresky_8k.hdr", skysphere_width, skysphere_height);
    //sycl::image<2> skysphere_hdr(skysphere_data.data(),
    //                             sycl::image_channel_order::rgba,
    //                             sycl::image_channel_type::fp32,
    //                             sycl::range<2>(skysphere_height, skysphere_width));

    std::cout << "[" << width << "x" << height << "]: " << RENDER_KERNEL_ITERATIONS * SAMPLES_PER_KERNEL << " samples" << std::endl << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    for (int kernel_iteration = 0; kernel_iteration < RENDER_KERNEL_ITERATIONS; kernel_iteration++)
    {
        try
        {
            queue.submit([&](sycl::handler& handler) {
                auto image_buffer_access = image_buffer.get_access<sycl::access::mode::write>(handler);
                auto triangle_buffer_access = triangle_buffer.get_access<sycl::access::mode::read>(handler);
                auto materials_buffer_access = materials_buffer.get_access<sycl::access::mode::read>(handler);
                auto emissive_triangle_indices_buffer_access = emissive_triangle_indices_buffer.get_access<sycl::access::mode::read>(handler);
                auto materials_indices_buffer_access = materials_indices_buffer.get_access<sycl::access::mode::read>(handler);
                //auto sphere_buffer_access = sphere_buffer.get_access<sycl::access::mode::read>(handler);
                auto bvh_nodes_access = bvh_nodes_buffer.get_access<sycl::access::mode::read>(handler);
                auto bvh_plane_normals = bvh_plane_normals_buffer.get_access<sycl::access::mode::read, sycl::access::target::device>(handler);
                /*auto skysphere_accessor = sycl::accessor<sycl::float4, 2, sycl::access::mode::read, sycl::access::target::image>(skysphere_hdr, handler);
                sycl::sampler skysphere_sampler(sycl::coordinate_normalization_mode::unnormalized, sycl::addressing_mode::clamp, sycl::filtering_mode::linear);*/

                const auto global_range = sycl::range<2>(width, height);
                const auto local_range = sycl::range<2>(TILE_SIZE_X, TILE_SIZE_Y);
                const auto coordinates_indices = sycl::nd_range<2>(global_range, local_range);

                sycl::stream debug_out_stream(2048, 256, handler);

                auto render_kernel = RenderKernel(width, height, kernel_iteration,
                    image_buffer_access,
                    triangle_buffer_access,
                    materials_buffer_access,
                    emissive_triangle_indices_buffer_access,
                    materials_indices_buffer_access,
                    //sphere_buffer_access,
                    bvh_nodes_access,
                    /*skysphere_accessor,
                    skysphere_sampler,*/
                    debug_out_stream);
                render_kernel.set_camera(Camera(45, Translation(0, 1, 3.5)));
                render_kernel.set_bvh_plane_normals(bvh_plane_normals);

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

//#include <sycl/sycl.hpp>
//#include <vector>
//#include <iostream>
//
//constexpr static size_t N{1024};
//
//int main() {
//  std::vector<int> v(N);
//  for (size_t i{0}; i<N; ++i) {
//    v[i] = i;
//  }
//
//  sycl::buffer buf{v};
//  try {
//      sycl::queue q{sycl::default_selector_v};
//      std::cout << "Running on device: "
//          << q.get_device().get_info<sycl::info::device::name>() << "\n";
//
//      q.submit([&](sycl::handler& cgh) {
//          auto acc{ buf.get_access(cgh,sycl::read_write) };
//          cgh.parallel_for(N, [=](sycl::id<1> id) {
//              int value{ acc[id] };
//              acc[id] = 2 * value;
//              acc[id] += 1;
//              });
//      });
//  }
//  catch (sycl::exception const& e) {
//      std::cout << "An exception is caught for vector add.\n";
//      std::cout << e.what() << std::endl;
//      std::terminate();
//  }
//
//  auto acc{buf.get_host_access(sycl::read_only)};
//  for (size_t i : {0, 63, 255, 511, 1023}) {
//    std::cout << "v[" << i << "] = " << acc[i] << std::endl;
//  }
//
//  return 0;
//}
