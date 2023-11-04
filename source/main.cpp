#include <iostream>
#include <chrono>
#include <cmath>

#include <OpenImageDenoise/oidn.hpp>
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

/*
 * A blend factor of 1 gives only the noisy image. 0 only the denoised image
 */
void oidn_denoise(Image& image, Image& output, float blend_factor)
{
    // Create an Open Image Denoise device
    oidn::DeviceRef device = oidn::newDevice(); // CPU or GPU if available
    device.commit();


    // Create buffers for input/output images accessible by both host (CPU) and device (CPU/GPU)
    int width = image.width();
    int height = image.height();

    oidn::BufferRef colorBuf = device.newBuffer(width * height * 3 * sizeof(float));
    // Create a filter for denoising a beauty (color) image using optional auxiliary images too
    // This can be an expensive operation, so try no to create a new filter for every image!
    oidn::FilterRef filter = device.newFilter("RT"); // generic ray tracing filter
    filter.setImage("color", colorBuf, oidn::Format::Float3, width, height); // beauty
    filter.setImage("output", colorBuf, oidn::Format::Float3, width, height); // denoised beauty
    filter.set("hdr", true); // beauty image is HDR
    filter.commit();
    // Fill the input image buffers
    float* colorPtr = (float*)colorBuf.getData();
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            int index = y * width + x;

            colorPtr[index * 3 + 0] = image[index].r;
            colorPtr[index * 3 + 1] = image[index].g;
            colorPtr[index * 3 + 2] = image[index].b;
        }
    // Filter the beauty image

    filter.execute();

    float* denoised_ptr = (float*)colorBuf.getData();
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            int index = y * width + x;

            Color color = blend_factor * Color(denoised_ptr[index * 3 + 0], denoised_ptr[index * 3 + 1], denoised_ptr[index * 3 + 2]) 
                + (1.0f - blend_factor) * image[index];
            color.a = 1.0f;

            output[index] = color;
        }

    const char* errorMessage;
    if (device.getError(errorMessage) != oidn::Error::None)
        std::cout << "Error: " << errorMessage << std::endl;
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

    const int width = 2000;
    const int height = 2000;

    Image image(width, height);

    ParsedOBJ parsed_obj = Utils::parse_obj("../SYCL-ray-tracing/data/OBJs/pbrt_dragon.obj");

    //Sphere sphere = add_sphere_to_scene(parsed_obj, Point(0.3275, 0.7, 0.3725), 0.2, SimpleMaterial {Color(0.0f), Color(1.0f, 0.71, 0.29), 1.0f, 0.4f}, parsed_obj.triangles.size());
    //std::vector<Sphere> spheres = { sphere };
    std::vector<Sphere> spheres;

    BVH bvh(&parsed_obj.triangles);
    //FlattenedBVH flat_bvh = bvh.flatten();

    std::vector<Triangle> triangle_buffer = parsed_obj.triangles;
    std::vector<SimpleMaterial> materials_buffer = parsed_obj.materials;
    std::vector<int> emissive_triangle_indices_buffer = parsed_obj.emissive_triangle_indices;
    std::vector<int> materials_indices_buffer = parsed_obj.material_indices;
    std::vector<Sphere> sphere_buffer = spheres;

    int skysphere_width, skysphere_height;
    Image skysphere_data = Utils::read_image_float("../SYCL-ray-tracing/data/Skyspheres/evening_road_01_puresky_8k.hdr", skysphere_width, skysphere_height);

    std::cout << "[" << width << "x" << height << "]: " << SAMPLES_PER_KERNEL << " samples" << std::endl << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    auto render_kernel = RenderKernel(width, height,
        image,
        triangle_buffer,
        materials_buffer,
        emissive_triangle_indices_buffer,
        materials_indices_buffer,
        sphere_buffer,
        bvh,
        skysphere_data);
    //render_kernel.set_camera(Camera::CORNELL_BOX_CAMERA);
    render_kernel.set_camera(Camera::PBRT_DRAGON_CAMERA);

    render_kernel.render();

    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << "ms" << std::endl;

    Image image_denoised_08(image.width(), image.height());
    oidn_denoise(image, image_denoised_08, 0.8);

    write_image_png(image_denoised_08, "../TP_RT_output_good_08_exp0.75.png");

    return 0;
}
