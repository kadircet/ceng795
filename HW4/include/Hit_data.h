#pragma once
#ifndef HIT_DATA_H_
#define HIT_DATA_H_
#include "Vector3.h"
class Shape;
struct Hit_data {
  float t;
  const Shape* shape;
  Vector3 normal;
  // Material_data material;
};
#endif
