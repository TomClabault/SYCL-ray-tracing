//#include <iostream>
//#include <sycl/sycl.hpp>
//#include <chrono>
//#include <cmath>

//#include <rapidobj/rapidobj.hpp>

//#include <stb_image_write.h>

//#include "bvh.h"
//#include "camera.h"
//#include "image_io.h"
//#include "render_kernel.h"
//#include "simple_material.h"
//#include "tests.h"
//#include "triangle.h"
//#include "utils.h"

//#include "xorshift.h"

//int dichotomie(std::vector<float> bins, float random)
//{
//    return 0;
//}

//int main(int argc, char* argv[])
//{
//    regression_tests();
//    std::cout << std::endl;

//    std::vector<float> bins {4, 5, 8, 2, 1};
//    std::vector<float> bins2 (bins.size());
//    std::vector<int> results(bins.size());
//    float sum = 0;
//    std::for_each(bins.begin(), bins.end(), [&sum] (float element) {
//        sum += element;
//    });

//    float cumul = 0.0f;
//    for (int i = 0; i < bins.size(); i++)
//    {
//        bins2[i] = bins[i] / sum + cumul;
//        cumul += bins[i] / sum;
//    }

//    xorshift32_generator generator(58);

//    for (int j = 0; j<  1000000; j++)
//    {
//        float random = generator();

//        results[dichotomie(bins2, random)]++;

//        /*
//        for (int i = 0; i < bins.size(); i++)
//        {
//            if (random < bins2[i])
//            {
//                results[i]++;

//                break;
//            }
//        }
//        */
//    }

//    const int width = 1280;
//    const int height = 720;

//    sycl::queue queue {sycl::gpu_selector_v};
//    std::cout << "Using " << queue.get_device().get_info<sycl::info::device::name>() << std::endl;

//    Image image(width, height);

//    ParsedOBJ parsed_obj = Utils::parse_obj("../SYCL-ray-tracing/data/OBJs/cornell.obj");

//    BVH bvh(&parsed_obj.triangles);
//    FlattenedBVH flat_bvh = bvh.flatten();

//    sycl::buffer<Color> image_buffer(image.color_data(), image.width() * image.height());
//    sycl::buffer<Triangle> triangle_buffer(parsed_obj.triangles.data(), parsed_obj.triangles.size());
//    sycl::buffer<SimpleMaterial> materials_buffer(parsed_obj.materials.data(), parsed_obj.materials.size());
//    sycl::buffer<int> emissive_triangle_indices_buffer(parsed_obj.emissive_triangle_indices.data(), parsed_obj.emissive_triangle_indices.size());
//    sycl::buffer<int> materials_indices_buffer(parsed_obj.material_indices.data(), parsed_obj.material_indices.size());
//    sycl::buffer<FlattenedBVH::FlattenedNode> bvh_nodes_buffer(flat_bvh.get_nodes().data(), flat_bvh.get_nodes().size());
//    sycl::buffer<Vector> bvh_plane_normals_buffer(BVH::BoundingVolume::PLANE_NORMALS, BVHConstants::PLANES_COUNT);

//    int skysphere_width, skysphere_height;
//    std::vector<sycl::float4> skysphere_data = Utils::read_image_float("../SYCL-ray-tracing/data/Skyspheres/evening_road_01_puresky_8k.hdr", skysphere_width, skysphere_height);
//    sycl::image<2> skysphere_hdr(skysphere_data.data(),
//                                 sycl::image_channel_order::rgba,
//                                 sycl::image_channel_type::fp32,
//                                 sycl::range<2>(skysphere_height, skysphere_width));

//    std::cout << "[" << width << "x" << height << "]: " << RENDER_KERNEL_ITERATIONS * SAMPLES_PER_KERNEL << " samples" << std::endl << std::endl;

//    auto start = std::chrono::high_resolution_clock::now();
//    for (int kernel_iteration = 0; kernel_iteration < RENDER_KERNEL_ITERATIONS; kernel_iteration++)
//    {
//        queue.submit([&] (sycl::handler& handler) {
//            auto image_buffer_access = image_buffer.get_access<sycl::access::mode::write>(handler);
//            auto triangle_buffer_access = triangle_buffer.get_access<sycl::access::mode::read>(handler);
//            auto materials_buffer_access = materials_buffer.get_access<sycl::access::mode::read>(handler);
//            auto emissive_triangle_indices_buffer_access = emissive_triangle_indices_buffer.get_access<sycl::access::mode::read>(handler);
//            auto materials_indices_buffer_access = materials_indices_buffer.get_access<sycl::access::mode::read>(handler);
//            auto bvh_nodes_access = bvh_nodes_buffer.get_access<sycl::access::mode::read>(handler);
//            auto bvh_plane_normals = bvh_plane_normals_buffer.get_access<sycl::access::mode::read, sycl::access::target::constant_buffer>(handler);
//            auto skysphere_accessor = sycl::accessor<sycl::float4, 2, sycl::access::mode::read, sycl::access::target::image>(skysphere_hdr, handler);
//            sycl::sampler skysphere_sampler(sycl::coordinate_normalization_mode::unnormalized, sycl::addressing_mode::clamp, sycl::filtering_mode::linear);

//            const auto global_range = sycl::range<2>(width, height);
//            const auto local_range = sycl::range<2>(TILE_SIZE_X, TILE_SIZE_Y);
//            const auto coordinates_indices = sycl::nd_range<2>(global_range, local_range);

//            sycl::stream debug_out_stream(1024, 256, handler);

//            auto render_kernel = RenderKernel(width, height, kernel_iteration,
//                                              image_buffer_access,
//                                              triangle_buffer_access,
//                                              materials_buffer_access,
//                                              emissive_triangle_indices_buffer_access,
//                                              materials_indices_buffer_access,
//                                              bvh_nodes_access,
//                                              skysphere_accessor,
//                                              skysphere_sampler,
//                                              debug_out_stream);
//            render_kernel.set_camera(Camera(45, Translation(0, 1, 3.5)));
//            render_kernel.set_bvh_plane_normals(bvh_plane_normals);

//            handler.parallel_for(coordinates_indices, render_kernel);
//        }).wait();

//        std::cout << (float)(kernel_iteration + 1) / RENDER_KERNEL_ITERATIONS * 100.0f << "%" << std::endl;
//    }

//    queue.wait();

//    auto stop = std::chrono::high_resolution_clock::now();
//    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << "ms" << std::endl;

//    image_buffer.get_access<sycl::access::mode::read>();


//    write_image_png(image, "../TP_RT_output.png");

//    return 0;
//}

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

//==============================================================
// Vector Add is the equivalent of a Hello, World! sample for data parallel
// programs. Building and running the sample verifies that your development
// environment is setup correctly and demonstrates the use of the core features
// of SYCL. This sample runs on both CPU and GPU (or FPGA). When run, it
// computes on both the CPU and offload device, then compares results. If the
// code executes on both CPU and offload device, the device name and a success
// message are displayed. And, your development environment is setup correctly!
//
// For comprehensive instructions regarding SYCL Programming, go to
// https://software.intel.com/en-us/oneapi-programming-guide and search based on
// relevant terms noted in the comments.
//
// SYCL material used in the code sample:
// •	A one dimensional array of data.
// •	A device queue, buffer, accessor, and kernel.
//==============================================================
// Copyright © Intel Corporation
//
// SPDX-License-Identifier: MIT
// =============================================================
#include <sycl/sycl.hpp>
#include <vector>
#include <iostream>
#include <string>
#if FPGA_HARDWARE || FPGA_EMULATOR || FPGA_SIMULATOR
#include <sycl/ext/intel/fpga_extensions.hpp>
#endif

using namespace sycl;

// num_repetitions: How many times to repeat the kernel invocation
size_t num_repetitions = 1;
// Vector type and data size for this example.
size_t vector_size = 10000;
typedef std::vector<int> IntVector;

// Create an exception handler for asynchronous SYCL exceptions
static auto exception_handler = [](sycl::exception_list e_list) {
    for (std::exception_ptr const& e : e_list) {
        try {
            std::rethrow_exception(e);
        }
        catch (std::exception const& e) {
#if _DEBUG
            std::cout << "Failure" << std::endl;
#endif
            std::terminate();
        }
    }
    };

//************************************
// Vector add in SYCL on device: returns sum in 4th parameter "sum_parallel".
//************************************
void VectorAdd(queue& q, const IntVector& a_vector, const IntVector& b_vector,
    IntVector& sum_parallel) {
    // Create the range object for the vectors managed by the buffer.
    range<1> num_items{ a_vector.size() };

    // Create buffers that hold the data shared between the host and the devices.
    // The buffer destructor is responsible to copy the data back to host when it
    // goes out of scope.
    buffer a_buf(a_vector);
    buffer b_buf(b_vector);
    buffer sum_buf(sum_parallel.data(), num_items);

    for (size_t i = 0; i < num_repetitions; i++) {

        // Submit a command group to the queue by a lambda function that contains the
        // data access permission and device computation (kernel).
        q.submit([&](handler& h) {
            // Create an accessor for each buffer with access permission: read, write or
            // read/write. The accessor is a mean to access the memory in the buffer.
            accessor a(a_buf, h, read_only);
            accessor b(b_buf, h, read_only);

            // The sum_accessor is used to store (with write permission) the sum data.
            accessor sum(sum_buf, h, write_only, no_init);

            // Use parallel_for to run vector addition in parallel on device. This
            // executes the kernel.
            //    1st parameter is the number of work items.
            //    2nd parameter is the kernel, a lambda that specifies what to do per
            //    work item. The parameter of the lambda is the work item id.
            // SYCL supports unnamed lambda kernel by default.
            h.parallel_for(num_items, [=](auto i) { sum[i] = a[i] + b[i]; });
            });
    };
    // Wait until compute tasks on GPU done
    q.wait();
}

//************************************
// Initialize the vector from 0 to vector_size - 1
//************************************
void InitializeVector(IntVector& a) {
    for (size_t i = 0; i < a.size(); i++) a.at(i) = i;
}

//************************************
// Demonstrate vector add both in sequential on CPU and in parallel on device.
//************************************
int main(int argc, char* argv[]) {
    // Change num_repetitions if it was passed as argument
    if (argc > 2) num_repetitions = std::stoi(argv[2]);
    // Change vector_size if it was passed as argument
    if (argc > 1) vector_size = std::stoi(argv[1]);
    // Create device selector for the device of your interest.
#if FPGA_EMULATOR
  // Intel extension: FPGA emulator selector on systems without FPGA card.
    auto selector = sycl::ext::intel::fpga_emulator_selector_v;
#elif FPGA_SIMULATOR
  // Intel extension: FPGA simulator selector on systems without FPGA card.
    auto selector = sycl::ext::intel::fpga_simulator_selector_v;
#elif FPGA_HARDWARE
  // Intel extension: FPGA selector on systems with FPGA card.
    auto selector = sycl::ext::intel::fpga_selector_v;
#else
  // The default device selector will select the most performant device.
    auto selector = default_selector_v;
#endif

    // Create vector objects with "vector_size" to store the input and output data.
    IntVector a, b, sum_sequential, sum_parallel;
    a.resize(vector_size);
    b.resize(vector_size);
    sum_sequential.resize(vector_size);
    sum_parallel.resize(vector_size);

    // Initialize input vectors with values from 0 to vector_size - 1
    InitializeVector(a);
    InitializeVector(b);

    try {
        queue q(selector, exception_handler);

        // Print out the device information used for the kernel code.
        std::cout << "Running on device: "
            << q.get_device().get_info<info::device::name>() << "\n";
        std::cout << "Vector size: " << a.size() << "\n";

        // Vector addition in SYCL
        VectorAdd(q, a, b, sum_parallel);
    }
    catch (exception const& e) {
        std::cout << "An exception is caught for vector add.\n";
        std::cout << e.what() << std::endl;
        std::terminate();
    }

    // Compute the sum of two vectors in sequential for validation.
    for (size_t i = 0; i < sum_sequential.size(); i++)
        sum_sequential.at(i) = a.at(i) + b.at(i);

    // Verify that the two vectors are equal.  
    for (size_t i = 0; i < sum_sequential.size(); i++) {
        if (sum_parallel.at(i) != sum_sequential.at(i)) {
            std::cout << "Vector add failed on device.\n";
            return -1;
        }
    }

    int indices[]{ 0, 1, 2, (static_cast<int>(a.size()) - 1) };
    constexpr size_t indices_size = sizeof(indices) / sizeof(int);

    // Print out the result of vector add.
    for (int i = 0; i < indices_size; i++) {
        int j = indices[i];
        if (i == indices_size - 1) std::cout << "...\n";
        std::cout << "[" << j << "]: " << a[j] << " + " << b[j] << " = "
            << sum_parallel[j] << "\n";
    }

    a.clear();
    b.clear();
    sum_sequential.clear();
    sum_parallel.clear();

    std::cout << "Vector add successfully completed on device.\n";
    return 0;
}