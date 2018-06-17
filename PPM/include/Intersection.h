#pragma once
#ifndef INTERSECTION_H_
#define INTERSECTION_H_
#include "Vector3.h"
class Shape;
class Intersection {
 public:
  Intersection()
      : t(std::numeric_limits<float>::infinity()),
        shape(nullptr),
        // is_light_object(false),
        normal(0.0f)
  // u(0.0f),
  // v(0.0f),
  // perlin_value(0.0f),
  // radiance(0.0f)
  {}

  float t;
  const Shape* shape;
  Vector3 normal;
  // float u, v;          // for textures
  // float perlin_value;  // for perlin noise
  // bool is_light_object;
  // For light objects
  // Vector3 radiance;
};
#endif
