#pragma once
#ifndef SCENE_H_
#define SCENE_H_
#include <string>
#include <vector>
#include "Camera.h"
#include "Material.h"
#include "Mesh.h"
#include "Photographic_tmo.h"
#include "Shape.h"
#include "Transformation.h"
#include "Vector3.h"
#include "Vertex.h"
class Scene {
 public:
  float shadow_ray_epsilon;
  Shape* bvh;
  std::vector<Mesh*> meshes;
  int max_recursion_depth;
  Vector3 ambient_light;
  std::vector<Scaling> scaling_transformations;
  std::vector<Translation> translation_transformations;
  std::vector<Rotation> rotation_transformations;
  std::vector<Camera> cameras;
  std::vector<Material> materials;
  std::vector<Vertex> vertex_data;
  Scene(const std::string& file_name);
  ~Scene();
  inline const Vertex& get_vertex_at(int index) const {
    return vertex_data[index];
  }
};
#endif
