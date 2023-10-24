#ifndef BVH_H
#define BVH_H

#include <array>
#include <cmath>
#include <limits>
#include <queue>

#include "bvh_constants.h"
#include "flattened_bvh.h"
#include "triangle.h"
#include "ray.h"

class FlattenedBVH;

class BVH
{
public:
    struct BoundingVolume
    {
        static Vector PLANE_NORMALS[BVHConstants::PLANES_COUNT];

        float _d_near[BVHConstants::PLANES_COUNT];
        float _d_far[BVHConstants::PLANES_COUNT];

        BoundingVolume()
        {
            for (int i = 0; i < BVHConstants::PLANES_COUNT; i++)
            {
                _d_near[i] = INFINITY;
                _d_far[i] = -INFINITY;
            }
        }

        static void triangle_volume(const Triangle& triangle, float d_near[BVHConstants::PLANES_COUNT], float d_far[BVHConstants::PLANES_COUNT])
        {
            for (int i = 0; i < BVHConstants::PLANES_COUNT; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    float dist = dot(BoundingVolume::PLANE_NORMALS[i], Vector(triangle[j]));

                    d_near[i] = std::min(d_near[i], dist);
                    d_far[i] = std::max(d_far[i], dist);
                }
            }
        }

        void extend_volume(const float d_near[BVHConstants::PLANES_COUNT], const float d_far[BVHConstants::PLANES_COUNT])
        {
            for (int i = 0; i < BVHConstants::PLANES_COUNT; i++)
            {
                _d_near[i] = std::min(_d_near[i], d_near[i]);
                _d_far[i] = std::max(_d_far[i], d_far[i]);
            }
        }

        void extend_volume(const BoundingVolume& volume)
        {
            extend_volume(volume._d_near, volume._d_far);
        }

        void extend_volume(const Triangle& triangle)
        {
            float d_near[BVHConstants::PLANES_COUNT];
            float d_far[BVHConstants::PLANES_COUNT];

            for (int i = 0; i < BVHConstants::PLANES_COUNT; i++)
            {
                d_near[i] = INFINITY;
                d_far[i] = -INFINITY;
            }

            triangle_volume(triangle, d_near, d_far);
            extend_volume(d_near, d_far);
        }

        static bool intersect(float d_near[BVHConstants::PLANES_COUNT],
                              float d_far[BVHConstants::PLANES_COUNT],
                              float denoms[BVHConstants::PLANES_COUNT],
                              float numers[BVHConstants::PLANES_COUNT])
        {
            float t_near = -INFINITY;
            float t_far = INFINITY;

            for (int i = 0; i < BVHConstants::PLANES_COUNT; i++)
            {
                float denom = denoms[i];
                if (denom == 0.0)
                    continue;

                //inverse denom to avoid division
                float d_near_i = (d_near[i] - numers[i]) / denom;
                float d_far_i = (d_far[i] - numers[i]) / denom;
                if (denom < 0)
                    std::swap(d_near_i, d_far_i);

                t_near = std::max(t_near, d_near_i);
                t_far = std::min(t_far, d_far_i);

                if (t_far < t_near)
                    return false;
            }

            return true;
        }

        /**
         * @params denoms Precomputed denominators
         */
        bool intersect(float& t_near, float& t_far, float* denoms, float* numers) const
        {
            t_near = -INFINITY;
            t_far = INFINITY;

            for (int i = 0; i < BVHConstants::PLANES_COUNT; i++)
            {
                float denom = denoms[i];
                if (denom == 0.0)
                    continue;

                //inverse denom to avoid division
                float d_near_i = (_d_near[i] - numers[i]) / denom;
                float d_far_i = (_d_far[i] - numers[i]) / denom;
                if (denom < 0)
                    std::swap(d_near_i, d_far_i);

                t_near = std::max(t_near, d_near_i);
                t_far = std::min(t_far, d_far_i);

                if (t_far < t_near)
                    return false;
            }

            return true;
        }
    };

    struct OctreeNode
    {
        struct QueueElement
        {
            QueueElement(const BVH::OctreeNode* node, float t_near) : _node(node), _t_near(t_near) {}

            bool operator > (const QueueElement& a) const
            {
                return _t_near > a._t_near;
            }

            const OctreeNode* _node;//Reference on the node

            float _t_near;//Intersection distance used to order the elements in the priority queue used
            //by the OctreeNode to compute the intersection with a ray
        };

        OctreeNode(Point min, Point max) : _min(min), _max(max) {}
        ~OctreeNode()
        {
            if (_is_leaf)
                return;
            else
            {
                for (int i = 0; i < 8; i++)
                    delete _children[i];
            }
        }

        /*
          * Once the objects have been inserted in the hierarchy, this function computes
          * the bounding volume of all the node in the hierarchy
          */
        BoundingVolume compute_volume()
        {
            if (_is_leaf)
                for (const Triangle* triangle : _triangles)
                    _bounding_volume.extend_volume(*triangle);
            else
                for (int i = 0; i < 8; i++)
                    _bounding_volume.extend_volume(_children[i]->compute_volume());

            return _bounding_volume;
        }

        void create_children(int max_depth, int leaf_max_obj_count)
        {
            float middle_x = (_min.x + _max.x) / 2;
            float middle_y = (_min.y + _max.y) / 2;
            float middle_z = (_min.z + _max.z) / 2;

            _children[0] = new OctreeNode(_min, Point(middle_x, middle_y, middle_z));
            _children[1] = new OctreeNode(Point(middle_x, _min.y, _min.z), Point(_max.x, middle_y, middle_z));
            _children[2] = new OctreeNode(_min + Point(0, middle_y, 0), Point(middle_x, _max.y, middle_z));
            _children[3] = new OctreeNode(Point(middle_x, middle_y, _min.z), Point(_max.x, _max.y, middle_z));
            _children[4] = new OctreeNode(_min + Point(0, 0, middle_z), Point(middle_x, middle_y, _max.z));
            _children[5] = new OctreeNode(Point(middle_x, _min.y, middle_z), Point(_max.x, middle_y, _max.z));
            _children[6] = new OctreeNode(_min + Point(0, middle_y, middle_z), Point(middle_x, _max.y, _max.z));
            _children[7] = new OctreeNode(Point(middle_x, middle_y, middle_z), Point(_max.x, _max.y, _max.z));
        }

        void insert(Triangle* triangle, int current_depth, int max_depth, int leaf_max_obj_count)
        {
            bool depth_exceeded = current_depth == max_depth;

            if (_is_leaf || depth_exceeded)
            {
                _triangles.push_back(triangle);

                if (_triangles.size() > leaf_max_obj_count && !depth_exceeded)
                {
                    _is_leaf = false;//This node isn't a leaf anymore

                    create_children(max_depth, leaf_max_obj_count);

                    for (Triangle* triangle : _triangles)
                        insert_to_children(triangle, current_depth, max_depth, leaf_max_obj_count);

                    _triangles.clear();
                    _triangles.shrink_to_fit();
                }
            }
            else
                insert_to_children(triangle, current_depth, max_depth, leaf_max_obj_count);

        }

        void insert_to_children(Triangle* triangle, int current_depth, int max_depth, int leaf_max_obj_count)
        {
            Point bbox_centroid = triangle->bbox_centroid();

            float middle_x = (_min.x + _max.x) / 2;
            float middle_y = (_min.y + _max.y) / 2;
            float middle_z = (_min.z + _max.z) / 2;

            int octant_index = 0;

            if (bbox_centroid.x > middle_x) octant_index += 1;
            if (bbox_centroid.y > middle_y) octant_index += 2;
            if (bbox_centroid.z > middle_z) octant_index += 4;

            _children[octant_index]->insert(triangle, current_depth + 1, max_depth, leaf_max_obj_count);
        }

        bool intersect(const Ray& ray, HitInfo& hit_info) const
        {
            float trash;

            float denoms[BVHConstants::PLANES_COUNT];
            float numers[BVHConstants::PLANES_COUNT];

            for (int i = 0; i < BVHConstants::PLANES_COUNT; i++)
            {
                denoms[i] = dot(BoundingVolume::PLANE_NORMALS[i], ray.direction);
                numers[i] = dot(BoundingVolume::PLANE_NORMALS[i], Vector(ray.origin));
            }

            return intersect(ray, hit_info, trash, denoms, numers);
        }

        bool intersect(const Ray& ray, HitInfo& hit_info, float& t_near, float* denoms, float* numers) const
        {
            float t_far, trash;

            if (!_bounding_volume.intersect(trash, t_far, denoms, numers))
                return false;

            if (_is_leaf)
            {
                for (Triangle* triangle : _triangles)
                {
                    HitInfo local_hit_info;
                    if (triangle->intersect(ray, local_hit_info))
                        if (local_hit_info.t < hit_info.t || hit_info.t == -1)
                            hit_info = local_hit_info;
                }

                t_near = hit_info.t;

                return t_near > 0;
            }

            std::priority_queue<QueueElement, std::vector<QueueElement>, std::greater<QueueElement>> intersection_queue;
            for (int i = 0; i < 8; i++)
            {
                float inter_distance;
                if (_children[i]->_bounding_volume.intersect(inter_distance, t_far, denoms, numers))
                    intersection_queue.emplace(QueueElement(_children[i], inter_distance));
            }

            bool intersection_found = false;
            float closest_inter = 100000000, inter_distance = 100000000;
            while (!intersection_queue.empty())
            {
                QueueElement top_element = intersection_queue.top();
                intersection_queue.pop();

                if (top_element._node->intersect(ray, hit_info, inter_distance, denoms, numers))
                {
                    closest_inter = std::min(closest_inter, inter_distance);
                    intersection_found = true;

                    //If we found an intersection that is closer than
                    //the next element in the queue, we can stop intersecting further
                    if (intersection_queue.empty() || closest_inter < intersection_queue.top()._t_near)
                    {
                        t_near = closest_inter;

                        return true;
                    }
                }
            }

            if (!intersection_found)
                return false;
            else
            {
                t_near = closest_inter;

                return true;
            }
        }

        std::vector<FlattenedBVH::FlattenedNode> flatten(int* current_node_index)
        {
//            std::vector<FlattenedBVH::FlattenedNode> nodes;

//            FlattenedBVH::FlattenedNode my_node;
//            for (int i = 0; i < BVHConstants::PLANES_COUNT; i++)
//            {
//                my_node.d_far[i] = _bounding_volume._d_far[i];
//                my_node.d_near[i] = _bounding_volume._d_near[i];
//            }

//            if (_is_leaf)
//            {
//                my_node.is_leaf = true;

//                for (int i = 0; i < 8; i++)
//                    my_node.children[i] = -1;

//                my_node.nb_triangles = _triangles.size();
//                for (int i = 0; i < _triangles.size(); i++)
//                    my_node.triangles_indices[i] = _triangles[i];

//                nodes.push_back(my_node);
//            }
//            else
//            {
//                for (int i = 0; i < 8; i++)
//                {
//                    (*current_node_index)++;
//                    std::vector<FlattenedBVH::FlattenedNode> child_nodes = _children[i]->flatten(current_node_index);
//                    my_node.children[i] = *current_node_index;

//                    nodes.insert(nodes.begin(), child_nodes.begin(), child_nodes.end());
//                }

//                nodes.push_back(my_node);
//            }

//            return nodes;

            return std::vector<FlattenedBVH::FlattenedNode>();
        }

        //If this node has been subdivided (and thus cannot accept any triangles),
        //this boolean will be set to false
        bool _is_leaf = true;

        std::vector<Triangle*> _triangles;
        BVH::OctreeNode* _children[8];

        Point _min, _max;
        BVH::BoundingVolume _bounding_volume;
    };

public:
    BVH();
    BVH(std::vector<Triangle>* triangles, int max_depth = 10, int leaf_max_obj_count = 8);
    ~BVH();

    void operator=(BVH&& bvh);

    bool intersect(const Ray& ray, HitInfo& hit_info) const;
    FlattenedBVH flatten() const;

private:
    void build_bvh(int max_depth, int leaf_max_obj_count, Point min, Point max, const BoundingVolume& volume);

public:
    OctreeNode* _root;

    std::vector<Triangle>* _triangles;
};

#endif
