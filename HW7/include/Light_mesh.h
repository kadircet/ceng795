#pragma once
#ifndef LIGHT_MESH_H_
#define LIGHT_MESH_H_
#include <vector>
#include "Light.h"
#include "Mesh.h"
#include "Mesh_triangle.h"
#include "Scene.h"
#include "Transformation.h"
#include "Vector3.h"
class Light_mesh : public Light, public Mesh {
 public:
  Light_mesh(const Scene* scene, int material_id,
             const Transformation& transformation,
             std::vector<Shape*>& triangles, const Vector3& radiance)
      : scene_(scene),
        radiance_(radiance),
        triangles_(triangles),
        Mesh(material_id, -1, triangles, transformation, 0.0f) {
    total_area_ = 0.0f;
    std::vector<float> pdf;
    const Matrix4x4& transformation_matrix =
        transformation.get_transformation_matrix();

    for (int i = 0; i < triangles.size(); i++) {
      Mesh_triangle* triangle = (Mesh_triangle*)triangles[i];
      Vector3 v_0 = transformation_matrix.multiply(
          scene_
              ->get_vertex_at(triangle->vertex_index_0 +
                              triangle->vertex_offset)
              .get_vertex_position());
      Vector3 v_1 = transformation_matrix.multiply(
          scene_
              ->get_vertex_at(triangle->vertex_index_1 +
                              triangle->vertex_offset)
              .get_vertex_position());
      Vector3 v_2 = transformation_matrix.multiply(
          scene_
              ->get_vertex_at(triangle->vertex_index_2 +
                              triangle->vertex_offset)
              .get_vertex_position());
      float area = (v_1 - v_0).cross(v_2 - v_0).length() / 2;
      total_area_ += area;
      pdf.push_back(area);
    }
    float current_cumulative_area = 0.0f;
    for (int i = 0; i < triangles.size(); i++) {
      current_cumulative_area += pdf[i];
      cdf_.push_back(current_cumulative_area / total_area_);
    }
  }
  bool intersect(const Ray& ray, Hit_data& hit_data,
                 bool culling) const override {
    const Matrix4x4& inverse_transformation_matrix =
        base_transform.get_inverse_transformation_matrix();
    Ray ray_local(inverse_transformation_matrix.multiply(ray.o),
                  inverse_transformation_matrix.multiply(ray.d, true),
                  ray.ray_type, ray.time);
    if (Mesh::intersect(ray_local, hit_data, culling)) {
      hit_data.normal = base_transform.get_normal_transformation_matrix()
                            .multiply(hit_data.normal, true)
                            .normalize();
      hit_data.is_light_object = true;
      hit_data.radiance = radiance_;
      return true;
    }
    return false;
  }
  Vector3 direction_and_distance(const Vector3& from_point,
                                 const Vector3& normal, float& distance,
                                 float& probability) const override;

  // Incoming radiance to the point from the light
  Vector3 incoming_radiance(const Vector3& from_point_to_light,
                            float probability) const override;

 private:
  Vector3 radiance_;
  // std::map<float, Shape*> triangles_;
  std::vector<Shape*> triangles_;
  std::vector<float> cdf_;
  float total_area_;
  const Scene* scene_;
};
#endif
