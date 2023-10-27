#ifndef PARSED_OBJ
#define PARSED_OBJ

#include "simple_material.h"
#include "triangle.h"

#include <vector>

struct ParsedOBJ
{
    std::vector<Triangle> triangles;
    std::vector<SimpleMaterial> materials;

    std::vector<int> emissive_triangle_indices;
    std::vector<int> material_indices;
};

#endif
