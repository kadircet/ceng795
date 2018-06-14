#pragma once
#ifndef LIGHT_SPHERE_H_
#define LIGHT_SPHERE_H_
#include "Light.h"
#include "Scene.h"
#include "Sphere.h"
#include "Transformation.h"
class Light_sphere : public Light, public Sphere {
 public:
  Light_sphere(const Scene* scene, const Vector3& center, float radius,
               int material_id, const Transformation& transformation,
               const Vector3& radiance)
      : radiance_(radiance),
        Sphere(scene, center, radius, material_id, -1, transformation,
               Vector3(0.0f)) {}
  bool intersect(const Ray& ray, Hit_data& hit_data) const override {
    if (Sphere::intersect(ray, hit_data)) {
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
};
#endif
