#pragma once
#ifndef SHAPE_H_
#define SHAPE_H_
#include "Hit_data.h"
class Ray;
class Shape {
 public:
  virtual Hit_data intersect(const Ray& ray) = 0;
};
#endif
