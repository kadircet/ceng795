#pragma once
#ifndef HIT_DATA_H_
#define HIT_DATA_H_
#include "Vector3.h"
class Shape;
struct Hit_data {
  float t;
  const Shape* shape;
  Vector3 normal;
  float u, v;          // for textures
  float perlin_value;  // for perlin noise
  // Material_data material;

  bool is_light_object;
  // For light objects
  Vector3 radiance;
};
#endif
