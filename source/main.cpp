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

int dichotomie(std::vector<float> bins, float random)
{
    return 0;
}

int main(int argc, char* argv[])
{
    regression_tests();
    std::cout << std::endl;

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

    const int width = 1280 / 4;
    const int height = 720 / 4;


    std::cout << "Reading OBJ..." << std::endl;
    //ParsedOBJ parsed_obj = Utils::parse_obj("../SYCL-ray-tracing/data/OBJs/pbrt_dragon.obj");
    ParsedOBJ parsed_obj = Utils::parse_obj("../SYCL-ray-tracing/data/OBJs/cornell_pbr.obj");

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
    Image skysphere_data = Utils::read_image_float("../SYCL-ray-tracing/data/Skyspheres/evening_road_01_puresky_8k.hdr", skysphere_width, skysphere_height);
    std::cout << "Importance Sampling Environment Map..." << std::endl;
    std::vector<ImageBin> skysphere_importance_bins;// = Utils::importance_split_skysphere(skysphere_data);

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
        skysphere_importance_bins);
    //render_kernel.set_camera(Camera::PBRT_DRAGON_CAMERA);
    render_kernel.set_camera(Camera::CORNELL_BOX_CAMERA);

    render_kernel.render();

    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << "ms" << std::endl;

    Image image_denoised_08 = Utils::OIDN_denoise(image_buffer, 0.8f);

    write_image_png(image_buffer, "../TP_RT_output_good_08_exp0.75.png");

    return 0;
}
