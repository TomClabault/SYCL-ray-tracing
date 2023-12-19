#include <iostream>
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

Sphere add_sphere_to_scene(ParsedOBJ& parsed_obj, const Point& center, float radius, const SimpleMaterial& material, int primitive_index)
{
    int material_index = parsed_obj.materials.size();

    parsed_obj.materials.push_back(material);
    parsed_obj.material_indices.push_back(material_index);

    Sphere sphere(center, radius, primitive_index);

    return sphere;
}

int main(int argc, char* argv[])
{
    regression_tests();
    std::cout << std::endl;

    const int width = 768;
    const int height = 768;


    std::cout << "Reading OBJ..." << std::endl;
    ParsedOBJ parsed_obj;
    //parsed_obj = Utils::parse_obj("../../data/cornell.obj");
    //parsed_obj = Utils::parse_obj("../SYCL-ray-tracing/data/OBJs/cornell_pbr.obj");
    //parsed_obj = Utils::parse_obj("../SYCL-ray-tracing/data/OBJs/ite-orb.obj");
    //parsed_obj = Utils::parse_obj("../SYCL-ray-tracing/data/OBJs/pbrt_dragon.obj");
    parsed_obj = Utils::parse_obj("../SYCL-ray-tracing/data/OBJs/ganesha_scene.obj");
    //parsed_obj = Utils::parse_obj("../SYCL-ray-tracing/data/OBJs/MIS.obj");

    //Sphere sphere = add_sphere_to_scene(parsed_obj, Point(0.3275, 0.7, 0.3725), 0.2, SimpleMaterial {Color(0.0f), Color(1.0f, 0.71, 0.29), 1.0f, 0.4f}, parsed_obj.triangles.size());
    //std::vector<Sphere> spheres = { sphere };
    std::vector<Sphere> spheres;

    BVH bvh(&parsed_obj.triangles);

    std::vector<Triangle> triangle_buffer = parsed_obj.triangles;
    std::vector<SimpleMaterial> materials_buffer = parsed_obj.materials;
    std::vector<int> emissive_triangle_indices_buffer = parsed_obj.emissive_triangle_indices;
    std::vector<int> materials_indices_buffer = parsed_obj.material_indices;
    std::vector<Sphere> sphere_buffer = spheres;

    int skysphere_width, skysphere_height;
    std::cout << "Reading Environment Map..." << std::endl;
    //Image skysphere_data = Utils::read_image_float("../SYCL-ray-tracing/data/Skyspheres/AllSkyFree_Sky_EpicGloriousPink_EquirectDebug.jpg", skysphere_width, skysphere_height);
    //Image skysphere_data = Utils::read_image_float("../SYCL-ray-tracing/data/Skyspheres/evening_road_01_puresky_8k.hdr", skysphere_width, skysphere_height);
    //Image skysphere_data = Utils::read_image_float("../SYCL-ray-tracing/data/Skyspheres/moonless_golf_8k.hdr", skysphere_width, skysphere_height);
    Image skysphere_data = Utils::read_image_float("../SYCL-ray-tracing/data/Skyspheres/rustig_koppie_puresky_8k.hdr", skysphere_width, skysphere_height);
    std::vector<float> env_map_cdf = Utils::compute_env_map_cdf(skysphere_data);

    std::cout << "[" << width << "x" << height << "]: " << SAMPLES_PER_KERNEL << " samples" << std::endl << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    Image image_buffer(width, height);
    auto render_kernel = RenderKernel(width, height,
        image_buffer,
        triangle_buffer,
        materials_buffer,
        emissive_triangle_indices_buffer,
        materials_indices_buffer,
        sphere_buffer,
        bvh,
        skysphere_data,
        env_map_cdf);
    //render_kernel.set_camera(Camera::CORNELL_BOX_CAMERA);
    render_kernel.set_camera(Camera::GANESHA_CAMERA);
    //render_kernel.set_camera(Camera::ITE_ORB_CAMERA);
    //render_kernel.set_camera(Camera::PBRT_DRAGON_CAMERA);
    //render_kernel.set_camera(Camera::MIS_CAMERA);

    render_kernel.render();

    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << "ms" << std::endl;

    Image image_denoised = Utils::OIDN_denoise(image_buffer, 1.0f);

    write_image_png(image_buffer, "../RT_output.png");
    write_image_png(image_denoised, "../RT_output_denoised.png");

    return 0;
}
