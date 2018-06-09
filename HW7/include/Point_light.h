#pragma once
#ifndef POINT_LIGHT_H
#define POINT_LIGHT_H
#include "Light.h"
class Point_light : public Light {
 public:
  Point_light(const Vector3& position, const Vector3& intensity);
  Vector3 direction_and_distance(const Vector3& from_point,
                                 float& distance) const override;
  Vector3 incoming_radiance(const Vector3& from_point_to_light) const override;

 private:
  Vector3 position_;
  Vector3 intensity_;
};
#endif
