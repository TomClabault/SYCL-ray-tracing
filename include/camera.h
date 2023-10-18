
#ifndef CAMERA_H
#define CAMERA_H

#include "mat.h"

struct Camera
{
    static inline const Transform DEFAULT_COORDINATES_SYSTEM = Transform(Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, -1), Vector(0, 0, 0));

    Camera() : view_matrix(DEFAULT_COORDINATES_SYSTEM)
    {
        fov = 45;
        full_fov_radians = fov / 180.0f * M_PI;
        fov_dist = 1.0f / std::tan(full_fov_radians / 2.0f);
    }

    /**
     * @brief Camera
     * @param full_fov In degrees
     * @param transformation
     */
    Camera(float full_fov, Transform transformation = Identity()) : view_matrix(transformation * DEFAULT_COORDINATES_SYSTEM)
    {
        fov = full_fov;
        full_fov_radians = fov / 180.0f * M_PI;
        fov_dist = 1.0f / std::tan(full_fov_radians / 2.0f);
    }

    Transform view_matrix;

    //Full FOV, not half
    float fov = 45;
    float full_fov_radians = fov / 180.0f * M_PI;
    float fov_dist = 1.0f / std::tan(full_fov_radians / 2.0f);
};

#endif
