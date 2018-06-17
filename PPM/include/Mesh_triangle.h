#pragma once
#ifndef MESH_TRIANGLE_H_
#define MESH_TRIANGLE_H_
#include <limits>
#include "Ray.h"
#include "Shape.h"
#include "Vector3.h"

class Scene;
enum Triangle_shading_mode { tsm_flat, tsm_smooth };

class Mesh_triangle : public Shape {
 public:
  int vertex_index_0, vertex_index_1, vertex_index_2;
  int vertex_offset;
  int texture_offset;
  Vector3 normal;
  int material_id;
  int texture_id;
  Triangle_shading_mode triangle_shading_mode;
  Mesh_triangle(const Scene* scene, int vertex_index_0, int vertex_index_1,
                int vertex_index_2, int vertex_offset, int texture_offset,
                int material_id, int texture_id, Triangle_shading_mode tsm);
  bool intersect(const Ray& ray, Intersection& intersection,
                 bool culling) const override;
  const Bounding_box& get_bounding_box() const override {
    return bounding_box_;
  }
  int get_material_id() const override { return material_id; }
  int get_texture_id() const override { return texture_id; }
  void print_debug(int indentation) const override {
    for (int index = 0; index < indentation; index++) {
      std::cout << "\t";
    }
    std::cout << "Mesh_triangle(" << vertex_index_0 << "," << vertex_index_1
              << "," << vertex_index_2 << "), material: " << material_id
              << "normal: " << normal << std::endl;
  }
  float get_surface_area() const;

 private:
  Bounding_box bounding_box_;
  const Scene* scene_;
  inline float determinant(const Vector3& col1, const Vector3& col2,
                           const Vector3& col3) const {
    return col1.x * (col2.y * col3.z - col3.y * col2.z) +
           col2.x * (col3.y * col1.z - col1.y * col3.z) +
           col3.x * (col1.y * col2.z - col2.y * col1.z);
  }
};
#endif
