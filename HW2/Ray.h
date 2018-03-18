#pragma once
#ifndef RAY_H_
#define RAY_H_
#include "Vector3.h"
class Ray {
 public:
  Vector3 o;
  Vector3 d;
  bool shadow;
  Ray(const Vector3& origin, const Vector3& direction, bool is_shadow = false)
      : o(origin), d(direction), shadow(is_shadow) {}
  inline Vector3 point_at(float t) const { return o + (t * d); }
};
#endif
