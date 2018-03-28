#pragma once
#ifndef SCENE_H_
#define SCENE_H_
#include <string>
#include <vector>
#include "Bounding_volume_hierarchy.h"
#include "Camera.h"
#include "Material.h"
#include "Mesh.h"
#include "Point_light.h"
#include "Shape.h"
#include "Sphere.h"
#include "Triangle.h"
#include "Vector3.h"
class Pixel;

class Scene {
 public:
  Vector3 background_color;
  float shadow_ray_epsilon;
  int max_recursion_depth;
  Vector3 ambient_light;
  Shape* bvh;
  std::vector<Camera> cameras;
  std::vector<Point_light> point_lights;
  std::vector<Material> materials;
  std::vector<Vector3> vertex_data;
  inline const Vector3& get_vertex_at(int index) const {
    return vertex_data[index];
  }

  Scene(const std::string& file_name);
  void render_image(int camera_index, Pixel* result, const int starting_row,
                    const int height_increase = 1) const;
  ~Scene();

 private:
  Vector3 trace_ray(const Ray& ray, int max_recursion_depth) const;
  bool refract_ray(const Vector3& direction_unit, const Vector3& normal,
                   const float refraction_index, Vector3& transmitted_d) const;
};
#endif
