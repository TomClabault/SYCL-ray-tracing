if (CMAKE_HOST_WIN32)
    # need at least CMake 3.23 for IntelLLVM support of IntelDPCPP package on Windows
    cmake_minimum_required(VERSION 3.23)
else()
    # CMake 3.20.5 is the minimum recommended for IntelLLVM on Linux
    cmake_minimum_required(VERSION 3.20.5)
endif()

project(SYCL_RayTracing LANGUAGES CXX)

find_package(IntelDPCPP REQUIRED)

add_executable(SYCL_RayTracing
    ${SYCL_HEADERS}

    source/main.cpp
    source/mat.cpp
    source/render_kernel.cpp

    include/camera.h
    include/color.h
    include/mat.h
    include/ray.h
    include/render_kernel.h
    include/triangle.h
    include/vec.h

    stbi/stb_image.h
    stbi/stb_image_write.h

    rapidobj/rapidobj.hpp)

set_property(TARGET SYCL_RayTracing PROPERTY CXX_STANDARD 17)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(SYCL_RayTracing PRIVATE Threads::Threads)

target_include_directories(SYCL_RayTracing PRIVATE "include/")
target_include_directories(SYCL_RayTracing PRIVATE "stbi/")
