#pragma once
#ifndef DIRECTIONAL_LIGHT_H_
#define DIRECTIONAL_LIGHT_H_
#include "Light.h"
class Directional_light : public Light {
 public:
  Directional_light(const Vector3& direction, const Vector3& radiance);
  Vector3 direction_and_distance(const Vector3& from_point, float& distance,
                                 float& probability) const override;
  Vector3 incoming_radiance(const Vector3& from_point_to_light,
                            float probability) const override;

 private:
  Vector3 direction_;
  Vector3 radiance_;
};
#endif
