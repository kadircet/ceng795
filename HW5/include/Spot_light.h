#pragma once
#ifndef SPOT_LIGHT_H_
#define SPOT_LIGHT_H_
#include "Light.h"
class Spot_light : public Light {
 public:
  Spot_light(const Vector3& position, const Vector3& intensity,
             const Vector3& direction, float coverage_angle_in_radians,
             float falloff_angle_in_radians);
  Vector3 direction_and_distance(const Vector3& from_point,
                                 float& distance) const override;
  Vector3 intensity(const Vector3& from_point_to_light) const override;

 private:
  Vector3 position_;
  Vector3 intensity_;
  Vector3 direction_;
  float cos_half_of_coverage_angle_;
  float cos_half_of_falloff_angle_;
};
#endif
