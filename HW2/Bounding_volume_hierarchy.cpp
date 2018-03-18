#include "Bounding_volume_hierarchy.h"
#include <iostream>
BVH::BVH(std::vector<Shape*>& objects, int start, int end, int dimension)
    : left(NULL), right(NULL) {
  for (int index = start; index < end; index++) {
    bounding_box.expand(objects[index]->get_bounding_box());
  }
  const float center = bounding_box.center[dimension];
  int mid_index = start;
  for (int index = start; index < end; index++) {
    const Shape* object = objects[index];
    if (object->get_bounding_box().center[dimension] < center) {
      std::swap(objects[index], objects[mid_index++]);
    }
  }
  if (mid_index == start || mid_index == end) {
    mid_index = start + ((end - start) / 2);
  }
  if (start + 1 == mid_index) {
    left = objects[start];
  } else {
    left = new BVH(objects, start, mid_index, (dimension + 1) % 3);
  }
  if (mid_index + 1 == end) {
    right = objects[mid_index];
  } else {
    right = new BVH(objects, mid_index, end, (dimension + 1) % 3);
  }
}

Hit_data BVH::intersect(const Ray& ray) const {
  Hit_data hit_data;
  hit_data.t = std::numeric_limits<float>::infinity();
  hit_data.shape = NULL;
  float bbox_t = bounding_box.intersect(ray);
  if (bbox_t < kEpsilon || bbox_t == kInf) {
    return hit_data;
  }
  Hit_data left_hit_data = left->intersect(ray);
  if (left_hit_data.t > kEpsilon && left_hit_data.t < hit_data.t) {
    hit_data = left_hit_data;
  }
  Hit_data right_hit_data = right->intersect(ray);
  if (right_hit_data.t > kEpsilon && right_hit_data.t < hit_data.t) {
    hit_data = right_hit_data;
  }

  return hit_data;
}
