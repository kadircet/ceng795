#pragma once
#ifndef MESH_H_
#define MESH_H_
#include <vector>
#include "Bounding_volume_hierarchy.h"
#include "Shape.h"
#include "Triangle.h"
class Mesh : public Shape {
 public:
  int material_id;
  int texture_id;
  Shape* bvh;
  bool intersect(const Ray& ray, Hit_data& hit_data) const override {
    return bvh->intersect(ray, hit_data);
  }
  int get_material_id() const override { return material_id; }
  const Bounding_box& get_bounding_box() const override {
    return bvh->get_bounding_box();
  }
  Mesh(int material_id, int texture_id, std::vector<Shape*>& triangles)
      : material_id(material_id), texture_id(texture_id), bvh(NULL) {
    bvh = BVH::create_bvh(triangles);
  }
  ~Mesh() {
    if (bvh) {
      delete bvh;
    }
  }
  void print_debug(int indentation) const override {
    for (int index = 0; index < indentation; index++) {
      std::cout << "\t";
    }
    std::cout << "Mesh->" << std::endl;
    bvh->print_debug(indentation);
  }
};
#endif
