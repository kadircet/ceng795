#pragma once
#ifndef SHAPE_H_
#define SHAPE_H_
#include "Hit_data.h"
class Ray;
class Shape {
 public:
  virtual Hit_data intersect(const Ray& ray) const = 0;
  virtual int get_material_id() const = 0;
};
#endif
