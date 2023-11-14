#include <iostream>
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

    std::vector<Color> image(width * height);
    std::vector<Triangle> triangles{ Triangle(Point(-1, 0, -2), Point(1, 0, -2), Point(0.5, 0.5, -2)) };

    std::cout << "[" << width << "x" << height << "]: " << SAMPLES_PER_KERNEL << " samples" << std::endl << std::endl;

    auto render_kernel = RenderKernel(width, height,
        image,
        triangles);
    render_kernel.set_camera(Camera(45, Translation(0, 1, 3.5)));
    render_kernel.render();

    write_image_png(image, width, height, "../TP_RT_output.png");

    return 0;
}
