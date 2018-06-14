#pragma once
#ifndef POINT_LIGHT_H
#define POINT_LIGHT_H
#include "Light.h"
class Point_light : public Light {
 public:
  Point_light(const Vector3& position, const Vector3& intensity);
  Vector3 direction_and_distance(const Vector3& from_point,
                                 const Vector3& normal, float& distance,
                                 float& probability) const override;
  Vector3 incoming_radiance(const Vector3& from_point_to_light,
                            float probability) const override;

 private:
  Vector3 position_;
  Vector3 intensity_;
};
#endif
