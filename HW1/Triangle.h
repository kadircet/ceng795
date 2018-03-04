#pragma once
#ifndef TRIANGLE_H_
#define TRIANGLE_H_
#include <limits>
#include "Ray.h"
#include "Shape.h"
#include "Vector3.h"

class Triangle : public Shape {
 public:
  Vector3 v_0, v_1, v_2;
  Vector3 normal;
  int material_id;
  Triangle(const Vector3& vertex_0, const Vector3& vertex_1,
           const Vector3& vertex_2, int material_id)
      : v_0(vertex_0), v_1(vertex_1), v_2(vertex_2), material_id(material_id) {
    normal = (v_1 - v_0).cross(v_2 - v_0).normalize();
  }
  Hit_data intersect(const Ray& ray) override {
    const Vector3 a_col1 = v_0 - v_1;
    const Vector3 a_col2 = v_0 - v_2;
    const Vector3& a_col3 = ray.d;
    const float det_a = determinant(a_col1, a_col2, a_col3);
    Hit_data hit_data;
    hit_data.t = std::numeric_limits<float>::infinity();
    hit_data.shape = NULL;
    if (det_a == 0) {
      return hit_data;
    }
    const Vector3& b = v_0 - ray.o;
    const float beta = determinant(b, a_col2, a_col3) / det_a;
    if (beta < -kEpsilon || beta > 1.0f + kEpsilon) return hit_data;
    const float gamma = determinant(a_col1, b, a_col3) / det_a;
    if (beta + gamma <= 1.0 + kEpsilon && beta >= -kEpsilon &&
        gamma >= -kEpsilon) {
      const float t = determinant(a_col1, a_col2, b) / det_a;
      hit_data.t = t;
      hit_data.shape = this;
      hit_data.normal = this->normal;
    }
    return hit_data;
  }

 private:
  inline float determinant(const Vector3& col1, const Vector3& col2,
                           const Vector3& col3) {
    return col1.x * (col2.y * col3.z - col3.y * col2.z) +
           col2.x * (col3.y * col1.z - col1.y * col3.z) +
           col3.x * (col1.y * col2.z - col2.y * col1.z);
  }
};
#endif
