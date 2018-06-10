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
             std::vector<Shape*>& triangles, const Vector3& radiance)
      : scene_(scene),
        radiance_(radiance),
        triangles_(triangles),
        Mesh(material_id, -1, triangles, Translation(0.0f, 0.0f, 0.0f), 0.0f) {
    total_area_ = 0.0f;
    std::vector<float> pdf;
    for (int i = 0; i < triangles.size(); i++) {
      Mesh_triangle* triangle = (Mesh_triangle*)triangles[i];
      float area = triangle->get_surface_area();
      total_area_ += area;
      pdf.push_back(area);
    }
    float current_cumulative_area = 0.0f;
    for (int i = 0; i < triangles.size(); i++) {
      current_cumulative_area += pdf[i];
      cdf_.push_back(current_cumulative_area / total_area_);
    }
  }
  bool intersect(const Ray& ray, Hit_data& hit_data) const override {
    if (Mesh::intersect(ray, hit_data)) {
      hit_data.is_light_object = true;
      hit_data.radiance = radiance_;
      return true;
    }
    return false;
  }
  Vector3 direction_and_distance(const Vector3& from_point, float& distance,
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
