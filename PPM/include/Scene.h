#pragma once
#ifndef SCENE_H_
#define SCENE_H_
#include <string>
#include <vector>
#include "Camera.h"
#include "Material.h"
#include "Photographic_tmo.h"
#include "Transformation.h"
#include "Vector3.h"
class Scene {
 public:
  float shadow_ray_epsilon;
  int max_recursion_depth;
  Vector3 ambient_light;
  std::vector<Scaling> scaling_transformations;
  std::vector<Translation> translation_transformations;
  std::vector<Rotation> rotation_transformations;
  std::vector<Camera> cameras;
  std::vector<Material> materials;
  Scene(const std::string& file_name);
  ~Scene();
};
#endif
