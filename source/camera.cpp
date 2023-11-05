#include "camera.h"

const Transform Camera::DEFAULT_COORDINATES_SYSTEM = Transform(Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, -1), Vector(0, 0, 0));
const Camera Camera::CORNELL_BOX_CAMERA = Camera(45, Translation(0, 1, 3.5));
const Camera Camera::PBRT_DRAGON_CAMERA = Camera(45, RotationX(-45) * Translation(0, -1, 10.5));