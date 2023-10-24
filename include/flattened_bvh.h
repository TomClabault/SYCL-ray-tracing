#ifndef FLATTENED_BVH_H
#define FLATTENED_BVH_H

#include "bvh_constants.h"
#include "hit_info.h"
#include "ray.h"
#include "triangle.h"

class FlattenedBVH
{
public:
    struct Stack
    {
        void push(int value) { stack[index++] = value; }
        int pop() { return stack[index--]; }
        bool empty() { return index == 0; }

        int stack[BVHConstants:: FLATTENED_BVH_MAX_STACK_SIZE];
        int index = 0;
    };

    struct FlattenedNode
    {
        bool intersect(float denoms[BVHConstants::PLANES_COUNT], float numers[BVHConstants::PLANES_COUNT]) const;

        //Extents of the planes of the bounding volume
        float d_near[BVHConstants::PLANES_COUNT];
        float d_far[BVHConstants::PLANES_COUNT];

        //Indices of the children in the m_nodes vector
        int children[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
        int triangles_indices[BVHConstants::MAX_TRIANGLES_PER_LEAF];

        int nb_triangles;
        int is_leaf;
    };

    bool intersect(const Ray& ray, HitInfo& hit_info, const std::vector<Triangle>& triangles) const;

    const std::vector<FlattenedNode>& get_nodes() const { return m_nodes; }
    std::vector<FlattenedNode>& get_nodes() { return m_nodes; }

private:
    std::vector<FlattenedNode> m_nodes;
};

#endif
