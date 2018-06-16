#pragma once
#ifndef SHAPE_H_
#define SHAPE_H_
#include "Bounding_box.h"
#include "Hit_data.h"
//#define CULLING_ENABLED
class Ray;
class Shape {
 public:
  virtual ~Shape() = 0;
  virtual const Bounding_box& get_bounding_box() const = 0;
  virtual bool intersect(const Ray& ray, Hit_data& hit_data,
                         bool culling) const = 0;
  virtual int get_material_id() const = 0;
  virtual int get_texture_id() const = 0;
  virtual void print_debug(int indent) const = 0;
};
#endif
