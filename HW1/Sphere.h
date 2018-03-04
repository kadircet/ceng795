#pragma once
#ifndef SPHERE_H_
#define SPHERE_H_
#include <cmath>
#include <limits>
#include "Ray.h"
#include "Shape.h"
#include "Vector3.h"

class Sphere : public Shape {
 public:
  Vector3 center;
  float radius;
  int material_id;

  Sphere(const Vector3& center, float radius, int material_id)
      : center(center), radius(radius), material_id(material_id) {}

  Hit_data intersect(const Ray& ray) override {
    Hit_data hit_data;
    hit_data.t = std::numeric_limits<float>::infinity();
    hit_data.shape = NULL;
    Vector3 center_to_origin = ray.o - center;
    const float a = ray.d.dot(ray.d);
    const float b = 2 * ray.d.dot(center_to_origin);
    const float c = center_to_origin.dot(center_to_origin) - radius * radius;
    const float determinant = b * b - 4 * a * c;
    if (determinant < -kEpsilon) {
      return hit_data;
    } else if (determinant < kEpsilon) {
      hit_data.t = -b / 2 * a;
      hit_data.normal = (ray.point_at(hit_data.t) - center).normalize();
      hit_data.shape = this;
    } else {
      const float t1 = (-b + sqrt(determinant)) / 2 * a;
      const float t2 = (-b - sqrt(determinant)) / 2 * a;
      // hit_data.t = fmin(t1, t2);
      if (t2 < .0) {
        hit_data.t = t1;
      } else {
        hit_data.t = t2;
      }
      hit_data.normal = (ray.point_at(hit_data.t) - center).normalize();
      hit_data.shape = this;
    }
    return hit_data;
  }
};
#endif
