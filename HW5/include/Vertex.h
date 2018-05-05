#pragma once
#ifndef VERTEX_H_
#define VERTEX_H_
#include "Vector3.h"

class Vertex {
 public:
  Vertex(const Vector3& position)
      : vertex_position_(position),
        vertex_normal_(0.0f),
        has_normal_(false),
        finalized_(false) {}
  inline const Vector3& get_vertex_position() const { return vertex_position_; }
  inline const Vector3& get_vertex_normal() const { return vertex_normal_; }
  inline bool has_vertex_normal() const { return has_normal_; }
  inline void add_vertex_normal(const Vector3& normal, float surface_area) {
    has_normal_ = true;
    vertex_normal_ += normal * surface_area;
  }
  inline void finalize_normal() {
    if (!finalized_) {
      vertex_normal_ = vertex_normal_.normalize();
      finalized_ = true;
    }
  }

 private:
  Vector3 vertex_position_;
  Vector3 vertex_normal_;
  bool has_normal_;
  bool finalized_;
};
#endif
