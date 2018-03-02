#pragma once
#ifndef CAMERA_H_
#define CAMERA_H_
#include <string>
#include "Image_plane"
#include "Vector3.h"
class Camera {
 public:
  Camera(const Vector3& up, const Vector3& gaze, const Vector3& position,
         const std::string& file, float left, float right, float bottom,
         float top, float distance, int image_width, int image_height)
      : e(position),
        file_name(file),
        image_plane(left, right, bottom, top, distance, image_width,
                    image_height) {
    w = -(gaze.normalize());
    u = up.normalize().cross(w).normalize();
    v = w.cross(u).normalize();
  }

 private:
  // implement a matrix frame for both basis and position?
  Vector3 u;  // side vector
  Vector3 v;  // up vector
  Vector3 w;  // forward vector
  Vector3 e;  // position
  Image_plane image_plane;
  const std::string file_name;
};
#endif
