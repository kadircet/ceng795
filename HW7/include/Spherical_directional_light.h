#pragma once
#ifndef SPHERICAL_DIRECTIONAL_LIGHT_H_
#define SPHERICAL_DIRECTIONAL_LIGHT_H_
#include <string>
#include <vector>
#include "Light.h"
#include "Vector3.h"

class Spherical_directional_light : public Light {
 public:
  Spherical_directional_light(const std::string& envmap_name);
  Vector3 direction_and_distance(const Vector3& from_point,
                                 const Vector3& normal, float& distance,
                                 float& probability) const override;
  Vector3 incoming_radiance(const Vector3& from_point_to_light,
                            float probability) const override;

 private:
  // only exr for now
  float* env_map_;
  int width_;
  int height_;
};
#endif
