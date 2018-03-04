#pragma once
#ifndef CAMERA_H_
#define CAMERA_H_
#include <string>
#include "Image_plane.h"
#include "Ray.h"
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

    top_left_corner = e - w * image_plane.distance + image_plane.left * u +
                      image_plane.top * v;
  }
  Ray calculate_ray_at(int x, int y) {
    // For now, origin is top_left, right handed camera
    const float s_u =
        (image_plane.right - image_plane.left) * (x + 0.5f) / image_plane.width;
    const float s_v = (image_plane.top - image_plane.bottom) * (y + 0.5f) /
                      image_plane.height;
    const Vector3& s = top_left_corner + s_u * u - s_v * v;
    return Ray(e, s - e);
  }

 private:
  // implement a matrix frame for both basis and position?
  Vector3 u;  // side vector
  Vector3 v;  // up vector
  Vector3 w;  // forward vector
  Vector3 e;  // position

  Vector3 top_left_corner;

  const std::string file_name;
  Image_plane image_plane;
};
#endif
