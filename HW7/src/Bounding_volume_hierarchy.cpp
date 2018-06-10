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

bool BVH::intersect(const Ray& ray, Hit_data& hit_data) const {
  float bbox_t = bounding_box.intersect(ray);
  if (bbox_t < 0.0f || bbox_t == kInf) {
    return false;
  }
  bool intersect = false;
  Hit_data left_hit_data;
  left_hit_data.t = std::numeric_limits<float>::infinity();
  left_hit_data.shape = NULL;
  left_hit_data.is_light_object = false;
  if (left->intersect(ray, left_hit_data) && left_hit_data.t > 0.0f &&
      left_hit_data.t < hit_data.t) {
    hit_data = left_hit_data;
    intersect = true;
  }
  Hit_data right_hit_data;
  right_hit_data.t = std::numeric_limits<float>::infinity();
  right_hit_data.shape = NULL;
  right_hit_data.is_light_object = false;
  if (right->intersect(ray, right_hit_data) && right_hit_data.t > 0.0f &&
      right_hit_data.t < hit_data.t) {
    intersect = true;
    hit_data = right_hit_data;
  }

  return intersect;
}
