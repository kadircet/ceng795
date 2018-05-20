#pragma once
#ifndef SCENE_H_
#define SCENE_H_
#include <string>
#include <vector>
#include "Area_light.h"
#include "BRDF.h"
#include "Bounding_volume_hierarchy.h"
#include "Camera.h"
#include "Directional_light.h"
#include "Material.h"
#include "Mesh.h"
#include "Mesh_triangle.h"
#include "Modified_phong_BRDF.h"
#include "Null_tmo.h"
#include "Phong_BRDF.h"
#include "Photographic_tmo.h"
#include "Point_light.h"
#include "Shape.h"
#include "Sphere.h"
#include "Spot_light.h"
#include "Texture.h"
#include "Transformation.h"
#include "Triangle.h"
#include "Vector3.h"
#include "Vertex.h"
class Pixel;

class Scene {
 public:
  Vector3 background_color;
  Texture* background_texture;
  float shadow_ray_epsilon;
  int max_recursion_depth;
  Vector3 ambient_light;
  Shape* bvh;
  std::vector<Mesh*> meshes;
  std::vector<Scaling> scaling_transformations;
  std::vector<Translation> translation_transformations;
  std::vector<Rotation> rotation_transformations;
  std::vector<Camera> cameras;
  std::vector<Light*> lights;
  std::vector<Material> materials;
  std::vector<Vertex> vertex_data;
  std::vector<Vector3> texture_coord_data;
  std::vector<Texture> textures;
  std::vector<BRDF*> brdfs;
  inline const Vertex& get_vertex_at(int index) const {
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
  void parse_ply_tinyply(std::string filename, std::vector<Vertex>& vertices,
                         std::vector<Shape*>& mesh_triangles, int vertex_offset,
                         int texture_offset, int material_id, int texture_id,
                         Triangle_shading_mode tsm) const;
};
#endif
