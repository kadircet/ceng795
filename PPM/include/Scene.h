#pragma once
#ifndef SCENE_H_
#define SCENE_H_
#include <mutex>
#include <string>
#include <thread>
#include <vector>
#include "Camera.h"
#include "Hit_point.h"
#include "Light.h"
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
  std::vector<Light*> lights;

  unsigned int num_hash;
  unsigned int num_photon;
  float hash_scale;
  std::vector<std::vector<Hit_point*>> hash_grid;
  std::vector<Hit_point*> hit_points;
  Bounding_box hit_point_bbox;

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
  void reset_hash_grid();
  void build_hash_grid(const int width, const int height);
  void eye_trace_lines(int index, int starting_row, int height_increase);
  void eye_trace(const Ray& ray, int depth, const Vector3& attenuation,
                 unsigned int pixel_index, float pixel_weight = 1.0f);
  void photon_trace(const Ray& ray, int depth, const Vector3& flux,
                    const Vector3& attenuation, int photon_id);
  void add_hit_point(Hit_point* hit_point) {
    std::lock_guard<std::mutex> lock(mutex_);
    hit_points.push_back(hit_point);
  };

 private:
  std::mutex mutex_;
  inline unsigned int hash(const int ix, const int iy, const int iz) {
    return (unsigned int)((ix * 73856093) ^ (iy * 19349663) ^ (iz * 83492791)) %
           num_hash;
  }
};
#endif
