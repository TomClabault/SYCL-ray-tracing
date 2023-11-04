//Already defined in image_io.cpp from gkit
//#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "color.h"
#include "image.h"
#include "rapidobj.hpp"
#include "utils.h"

#include <iostream>
#include <string>

ParsedOBJ Utils::parse_obj(const std::string& filepath)
{
    ParsedOBJ parsed_obj;

    rapidobj::Result rapidobj_result = rapidobj::ParseFile(filepath, rapidobj::MaterialLibrary::Default());
    if (rapidobj_result.error)
    {
        std::cout << "There was an error loading the OBJ file: " << rapidobj_result.error.code.message() << std::endl;
        std::cin.get();

        return parsed_obj;
    }

    rapidobj::Triangulate(rapidobj_result);

    const rapidobj::Array<float>& positions = rapidobj_result.attributes.positions;
    std::vector<Triangle>& triangle_buffer = parsed_obj.triangles;
    std::vector<int>& emissive_triangle_indices_buffer = parsed_obj.emissive_triangle_indices;
    std::vector<int>& materials_indices_buffer = parsed_obj.material_indices;
    for (rapidobj::Shape& shape : rapidobj_result.shapes)
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
            triangle_buffer.push_back(triangle);

            int mesh_triangle_index = i / 3;
            int triangle_material_index = mesh.material_ids[mesh_triangle_index];
            //Always go +1 because 0 is the default material. This also has for effect
            //to transform the -1 material index (on material) to 0 which is the default
            //material
            materials_indices_buffer.push_back(triangle_material_index + 1);

            rapidobj::Float3 emission;
            if (triangle_material_index == -1)
                emission = rapidobj::Float3{ 0, 0, 0 };
            else
                emission = rapidobj_result.materials[triangle_material_index].emission;
            if (emission[0] > 0 || emission[1] > 0 || emission[2] > 0)
            {
                //This is an emissive triangle
                //Using the buffer of all the triangles to get the global id of the emissive
                //triangle
                emissive_triangle_indices_buffer.push_back(triangle_buffer.size() - 1);
            }
        }
    }

    //Computing SimpleMaterials
    std::vector<SimpleMaterial>& materials_buffer = parsed_obj.materials;
    materials_buffer.push_back(SimpleMaterial{ Color(1.0f, 0.0f, 1.0f), Color(), 0.0f, 1.0f });
    for (const rapidobj::Material& material : rapidobj_result.materials)
        materials_buffer.push_back(SimpleMaterial {Color(material.emission), Color(material.diffuse), material.metallic, material.roughness});

    return parsed_obj;
}

Image Utils::read_image_float(const std::string& filepath, int& image_width, int& image_height)
{
    int channels;
    float* pixels = stbi_loadf(filepath.c_str(), &image_width, &image_height, &channels, 0);

    if(!pixels)
    {
        std::cout << "Error reading image " << filepath << std::endl;
        return Image();
    }

    Image output(image_width, image_height);
    //Using the pixels pointer as an iterator
    for (int y = 0; y < image_height; y++)
    {
        for (int x = 0; x < image_width; x++)
        {
            int index = y * image_width + x;
            output[index] = pow(Color(pixels[index * 3 + 0], pixels[index * 3 + 1], pixels[index * 3 + 2], 0.0f), 2.2f);
        }
    }

    return output;
}
