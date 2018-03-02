#pragma once
#ifndef SCENE_H_
#define SCENE_H_
#include <string>
#include <vector>
#include "Camera.h"
#include "Vector3i.cpp"
class Scene {
 public:
  Vector3i background_color;
  float shadow_ray_epsilon;
  int max_recursion_depth;
  Vector3 ambient_light;

  std::vector<Camera> cameras;
  Scene(const std::string& file_name);
};
#endif
