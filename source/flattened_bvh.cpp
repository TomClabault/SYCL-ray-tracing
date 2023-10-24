#include "flattened_bvh.h"

#include "bvh.h"

bool FlattenedBVH::FlattenedNode::intersect(float denoms[], float numers[]) const
{
    //TODO
    //return BVH::BoundingVolume::intersect(d_near, d_far, denoms, numers);
}

bool FlattenedBVH::intersect(const Ray& ray, HitInfo& hit_info, const std::vector<Triangle>& triangles) const
{
    hit_info.t = -1;

    FlattenedBVH::Stack stack;
    stack.push(0);//Pushing the root of the BVH

    float denoms[BVHConstants::PLANES_COUNT];
    float numers[BVHConstants::PLANES_COUNT];

    for (int i = 0; i < BVHConstants::PLANES_COUNT; i++)
    {
        denoms[i] = dot(BVH::BoundingVolume::PLANE_NORMALS[i], ray.direction);
        numers[i] = dot(BVH::BoundingVolume::PLANE_NORMALS[i], Vector(ray.origin));
    }

    while (!stack.empty())
    {
        int node_index = stack.pop();
        const FlattenedNode& node = m_nodes[node_index];

        if (node.intersect(denoms, numers))
        {
            if (node.is_leaf)
            {
                float closest_intersection_distance = std::numeric_limits<float>::max();
                for (int i = 0; i < node.nb_triangles; i++)
                {
                    int triangle_index = node.triangles_indices[i];

                    HitInfo local_hit_info;
                    if (triangles[triangle_index].intersect(ray, local_hit_info))
                    {
                        if (closest_intersection_distance > local_hit_info.t)
                        {
                            closest_intersection_distance = local_hit_info.t;
                            hit_info = local_hit_info;
                        }
                    }
                }
            }
            else
                for (int i = 0; i < 8; i++)
                    stack.push(node.children[i]);
        }
    }

    return hit_info.t > -1;
}
