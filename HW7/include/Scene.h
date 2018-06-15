#pragma once
#ifndef SCENE_H_
#define SCENE_H_
#include <string>
#include <vector>
#include "Area_light.h"
#include "BRDF.h"
#include "Blinn_phong_BRDF.h"
#include "Bounding_volume_hierarchy.h"
#include "Camera.h"
#include "Directional_light.h"
#include "Material.h"
#include "Mesh.h"
#include "Mesh_triangle.h"
#include "Modified_blinn_phong_BRDF.h"
#include "Modified_phong_BRDF.h"
#include "Phong_BRDF.h"
#include "Photographic_tmo.h"
#include "Point_light.h"
#include "Shape.h"
#include "Sphere.h"
#include "Spherical_directional_light.h"
#include "Spot_light.h"
#include "Texture.h"
#include "Torrance_sparrow_BRDF.h"
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
  Spherical_directional_light* spherical_directional_light;
  inline const Vertex& get_vertex_at(int index) const {
    return vertex_data[index];
  }

  Scene(const std::string& file_name);
  void render_image(int camera_index, Pixel* result, const int starting_row,
                    const int height_increase = 1) const;
  ~Scene();

 private:
  Vector3 send_ray(const Ray& ray, int recursion_level) const;
  Vector3 trace_ray(const Ray& ray, const Hit_data& hit_data,
                    int recursion_level) const;
  Vector3 trace_path(const Ray& ray, const Hit_data& hit_data,
                     int recursion_level) const;
  Vector3 reflect_ray(const Ray& ray, const Hit_data& hit_data,
                      int recursion_level) const;
  Vector3 refract_ray(const Ray& ray, const Hit_data& hit_data,
                      int recursion_level) const;
  Vector3 calculate_diffuse_and_specular_radiance(
      const Ray& ray, const Hit_data& hit_data,
      const Vector3& diffuse_constant) const;
  bool calculate_diffuse_constant(const Hit_data& hit_data,
                                  Vector3& diffuse_constant_out) const;
  bool calculate_transmission(const Vector3& direction_unit,
                              const Vector3& normal,
                              const float refraction_index,
                              Vector3& transmitted_d) const;
  void parse_ply_tinyply(std::string filename, std::vector<Vertex>& vertices,
                         std::vector<Shape*>& mesh_triangles, int vertex_offset,
                         int texture_offset, int material_id, int texture_id,
                         Triangle_shading_mode tsm);
};
#endif
