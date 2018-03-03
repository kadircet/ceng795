#pragma once
#ifndef SCENE_H_
#define SCENE_H_
#include <string>
#include <vector>
#include "Camera.h"
#include "Material.h"
#include "Point_light.h"
#include "Vector3i.cpp"

class Scene {
 public:
  Vector3i background_color;
  float shadow_ray_epsilon;
  int max_recursion_depth;
  Vector3 ambient_light;

  std::vector<Camera> cameras;
  std::vector<Point_light> point_lights;
  std::vector<Material> materials;
  std::vector<Vector3> vertex_data;
  Scene(const std::string& file_name);
};
#endif
