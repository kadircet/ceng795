#pragma once
#ifndef AREA_LIGHT_H
#define AREA_LIGHT_LIGHT_H
#include "Light.h"
class Area_light : public Light {
 public:
  Area_light(const Vector3& position, const Vector3& intensity,
             const Vector3& edge_vector_1, const Vector3& edge_vector_2);
  Vector3 direction_and_distance(const Vector3& from_point,
                                 float& distance) const override;
  Vector3 incoming_radiance(const Vector3& from_point_to_light) const override;

 private:
  Vector3 position_;
  Vector3 intensity_;
  Vector3 edge_vector_1_;
  Vector3 edge_vector_2_;
  Vector3 normal_;
};
#endif
