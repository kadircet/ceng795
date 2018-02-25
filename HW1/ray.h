#pragma once
#ifndef RAY_H_
#define RAY_H_
#include "Vector3.h"
class Ray {
 public:
  Vector3 o;
  Vector3 d;
  Ray(const Vector3& origin, const Vector3& direction)
      : o(origin), d(direction) {}
  inline Vector3 point_at(float t) const { return o + (t * d); }
};
#endif
