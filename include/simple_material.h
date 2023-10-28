#ifndef SIMPLE_MATERIAL_H
#define SIMPLE_MATERIAL_H

#include "color.h"

struct SimpleMaterial
{
    Color emission;
    Color diffuse;

    float metalness;
    float roughness;
};

#endif
