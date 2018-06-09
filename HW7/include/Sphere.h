#pragma once
#ifndef SPHERE_H_
#define SPHERE_H_
#include <cmath>
#include <limits>
#include "Shape.h"
#include "Transformation.h"
#include "Vector3.h"
class Ray;
class Scene;

class Sphere : public Shape {
 public:
  Vector3 center;
  float radius;
  int material_id;
  int texture_id;
  Vector3 velocity;

  Sphere(const Scene* scene, const Vector3& center, float radius, int material_id, int texture_id,
         const Transformation& transformation, const Vector3& velocity);
  
  int get_material_id() const override { return material_id; }
  int get_texture_id() const override { return texture_id; }
  const Bounding_box& get_bounding_box() const override {
    return bounding_box_;
  }
  bool intersect(const Ray& ray, Hit_data& hit_data) const override;
  void print_debug(int indentation) const override {
    for (int index = 0; index < indentation; index++) {
      std::cout << "\t";
    }
    std::cout << "Sphere: " << center << "," << radius
              << " material: " << material_id << std::endl;
  }

 private:
  Transformation transformation_;
  bool is_identity_;
  Bounding_box bounding_box_;
  const Scene* scene_;
  void get_uv(const Vector3& local_coordinates, float& u, float& v) const;
};
#endif
