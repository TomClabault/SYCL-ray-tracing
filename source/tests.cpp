#include "rapidobj.hpp"

#include "bvh.h"
#include "bvh_tests.h"
#include "flattened_bvh.h"
#include "ray.h"
#include "tests.h"
#include "triangle.h"

bool compare_points(const Point& a, const Point& b)
{
    return std::abs(a.x - b.x) < 1.0e-5f
            && std::abs(a.y - b.y) < 1.0e-5f && std::abs(a.z - b.z) < 1.0e-5f;
}

void test_bvh(BVH& bvh)
{
    std::cout << "[TESTS]: BVH MISSING RAYS --- ";
    for (const Ray& ray : bvh_test_rays_no_inter)
    {
        HitInfo closest_hit_info;
        if (bvh.intersect(ray, closest_hit_info))
        {
            std::cout << "Intersection found with the ray: " << ray << " but no intersection should have been found" << std::endl;
            std::cout << "Intersection t/point: " << closest_hit_info.t << "/" << ray.origin + closest_hit_info.t * ray.direction << std::endl;

            std::exit(-1);
        }
    }
    std::cout << "OK" << std::endl;

    std::cout << "[TESTS]: BVH HITTING RAYS --- ";
    for (int index = 0; index < bvh_test_rays_inter.size(); index++)
    {
        const Ray& ray = bvh_test_rays_inter[index];
        const Point& expected_point = bvh_test_rays_inter_result_points[index];

        HitInfo closest_hit_info;
        if (!bvh.intersect(ray, closest_hit_info))
        {
            std::cout << "No intersection was found with the ray: " << ray << " but an intersection should have been found" << std::endl;
            std::cout << "Expected intersection: " << expected_point << std::endl;

            std::exit(-1);
        }
        else
        {
            Point inter_point = ray.origin + closest_hit_info.t * ray.direction;
            if (!compare_points(inter_point, expected_point))
            {
                std::cout << "An intersection was found but the intersection is incorrect. Found point: " << inter_point << " | expected point: " << expected_point << std::endl;

                std::exit(-1);
            }
        }
    }
    std::cout << "OK" << std::endl;
}

void test_flattened_bvh(BVH& bvh)
{
    FlattenedBVH flat_bvh = bvh.flatten();

    std::cout << "[TESTS]: FLATTENED BVH MISSING RAYS --- ";
    for (const Ray& ray : bvh_test_rays_no_inter)
    {
        HitInfo closest_hit_info;
        if (flat_bvh.intersect(ray, closest_hit_info, *bvh._triangles))
        {
            std::cout << "Intersection found with the ray: " << ray << " but no intersection should have been found" << std::endl;
            std::cout << "Intersection t/point: " << closest_hit_info.t << "/" << ray.origin + closest_hit_info.t * ray.direction << std::endl;

            std::exit(-1);
        }
    }
    std::cout << "OK" << std::endl;

    std::cout << "[TESTS]: FLATTENED BVH HITTING RAYS --- ";
    for (int index = 0; index < bvh_test_rays_inter.size(); index++)
    {
        const Ray& ray = bvh_test_rays_inter[index];
        const Point& expected_point = bvh_test_rays_inter_result_points[index];

        HitInfo closest_hit_info;
        if (!flat_bvh.intersect(ray, closest_hit_info, *bvh._triangles))
        {
            std::cout << "No intersection was found with the ray: " << ray << " but an intersection should have been found" << std::endl;
            std::cout << "Expected intersection: " << expected_point << std::endl;

            std::exit(-1);
        }
        else
        {
            Point inter_point = ray.origin + closest_hit_info.t * ray.direction;
            if (!compare_points(inter_point, expected_point))
            {
                std::cout << "An intersection was found but the intersection is incorrect. Found point: " << inter_point << " | expected point: " << expected_point << std::endl;

                std::exit(-1);
            }
        }
    }
    std::cout << "OK" << std::endl;
}

void regression_tests()
{
    rapidobj::Result parsed_obj = rapidobj::ParseFile("../SYCL-ray-tracing/data/cornell.obj", rapidobj::MaterialLibrary::Default());
    if (parsed_obj.error)
    {
        std::cout << "There was an error loading the OBJ file: " << parsed_obj.error.code.message() << std::endl;

        std::exit(-1);
    }
    rapidobj::Triangulate(parsed_obj);

    const rapidobj::Array<float>& positions = parsed_obj.attributes.positions;
    std::vector<Triangle> triangle_host_buffer;
    for (rapidobj::Shape& shape : parsed_obj.shapes)
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
            triangle_host_buffer.push_back(triangle);
        }
    }

    BVH bvh(&triangle_host_buffer);
    test_bvh(bvh);
    test_flattened_bvh(bvh);
}
