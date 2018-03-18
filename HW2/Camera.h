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
         int number_of_samples, const std::string& filename, float left,
         float right, float bottom, float top, float distance, int image_width,
         int image_height)
      : e(position),
        number_of_samples_(number_of_samples),
        filename_(filename),
        image_plane_(left, right, bottom, top, distance, image_width,
                     image_height) {
    w = -(gaze.normalize());
    u = up.normalize().cross(w).normalize();
    v = w.cross(u).normalize();

    top_left_corner = e - w * image_plane_.distance + image_plane_.left * u +
                      image_plane_.top * v;
    s_u_constant =
        ((image_plane_.right - image_plane_.left) / image_plane_.width) * u;
    s_v_constant =
        ((image_plane_.top - image_plane_.bottom) / image_plane_.height) * v;
  }
  Ray calculate_ray_at(float x, float y) const {
    // For now, origin is top_left, right handed camera
    const Vector3& s =
        top_left_corner + (x + 0.5) * s_u_constant - (y + 0.5) * s_v_constant;
    return Ray(e, (s - e));
  }
  const Image_plane& get_image_plane() const { return image_plane_; }
  const std::string& get_filename() const { return filename_; }
  int get_number_of_samples() const { return number_of_samples_; }

 private:
  // implement a matrix frame for both basis and position?
  Vector3 u;  // side vector
  Vector3 v;  // up vector
  Vector3 w;  // forward vector
  Vector3 e;  // position

  Vector3 s_u_constant;
  Vector3 s_v_constant;
  Vector3 top_left_corner;

  int number_of_samples_;
  const std::string filename_;
  Image_plane image_plane_;
};
#endif
