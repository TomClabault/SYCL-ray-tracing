#ifndef UTILS_H
#define UTILS_H

#include "parsed_obj.h"

#include <string>
#include <sycl/sycl.hpp>

class Utils
{
public:
    static ParsedOBJ parse_obj(const std::string& filepath);
    static std::vector<sycl::float4> read_image_float(const std::string& filepath, int& image_width, int& image_height);
};

#endif
