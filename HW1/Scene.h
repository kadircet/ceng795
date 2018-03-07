#pragma once
#ifndef SCENE_H_
#define SCENE_H_
#include <string>
#include <vector>
#include "Camera.h"
#include "Material.h"
#include "Mesh.h"
#include "Point_light.h"
#include "Sphere.h"
#include "Triangle.h"
#include "Vector3.h"
#include "Vector3i.h"

class Scene {
 public:
  Vector3 background_color;
  float shadow_ray_epsilon;
  int max_recursion_depth;
  Vector3 ambient_light;

  std::vector<Camera> cameras;
  std::vector<Point_light> point_lights;
  std::vector<Material> materials;
  std::vector<Vector3> vertex_data;
  inline const Vector3& get_vertex_at(int index) const {
    return vertex_data[index];
  }

  std::vector<Mesh> meshes;
  std::vector<Triangle> triangles;
  std::vector<Sphere> spheres;
  Scene(const std::string& file_name);
  void render_image(int camera_index, Vector3i* result, const int starting_row,
                    const int height_increase = 0);

 private:
  Vector3 trace_ray(const Ray& ray, int max_recursion_depth);
};
#endif
