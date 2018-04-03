#pragma once
#ifndef RAY_H_
#define RAY_H_
#include "Vector3.h"
enum Ray_type {
  r_primary,
  r_shadow,
  r_reflection,
  r_refraction
};
class Ray {
 public:
  Vector3 o;
  Vector3 d;
  bool in_medium;
  Ray_type ray_type;
  Ray(const Vector3& origin, const Vector3& direction, Ray_type ray_type)
      : o(origin), d(direction), in_medium(false), ray_type(ray_type) {}
  inline Vector3 point_at(float t) const { return o + (t * d); }
};
#endif
