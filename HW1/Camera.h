#pragma once
#ifndef CAMERA_H_
#define CAMERA_H_
#include "Vector3.h"
class Camera {
 public:
  // implement a matrix frame for both basis and position?
  Vector3 u;  // side vector
  Vector3 v;  // up vector
  Vector3 w;  // forward vector
  Vector3 e;  // position
  Camera(const Vector3& up, const Vector3& gaze, const Vector3& position)
      : e(position) {
    w = -(gaze.normalize());
    u = up.normalize().cross(w).normalize();
    v = w.cross(u).normalize();
  }
};
#endif
