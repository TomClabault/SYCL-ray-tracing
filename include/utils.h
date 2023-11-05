#ifndef UTILS_H
#define UTILS_H

#include "parsed_obj.h"

#include <string>

class Utils
{
public:
    static ParsedOBJ parse_obj(const std::string& filepath);
    static Image read_image_float(const std::string& filepath, int& image_width, int& image_height);

    /*
     * A blend factor of 1 gives only the noisy image. 0 only the denoised image
     */
    static Image OIDN_denoise(const Image& noisy_image, float blend_factor);
};

#endif
