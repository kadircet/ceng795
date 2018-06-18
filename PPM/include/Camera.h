#pragma once
#ifndef CAMERA_H_
#define CAMERA_H_
#include <string>
#include "Ray.h"
#include "Tonemapping_operator.h"
#include "Vector3.h"
#include "tinyxml2.h"

struct Image_plane {
  float left, right, bottom, top;
  float distance;
  int width, height;
  Image_plane(float l, float r, float b, float t, float d, int image_width,
              int image_height)
      : left(l),
        right(r),
        bottom(b),
        top(t),
        distance(d),
        width(image_width),
        height(image_height) {}
};

class Camera {
 public:
  static void load_cameras_from_xml(tinyxml2::XMLElement* element,
                                    std::vector<Camera>& cameras);
  Camera(const Vector3& up, const Vector3& gaze, const Vector3& position,
         int number_of_samples, const std::string& filename, float left,
         float right, float bottom, float top, float distance, int image_width,
         int image_height, Tonemapping_operator* tmo, bool left_handed)
      : e(position),
        number_of_samples_(number_of_samples),
        filename_(filename),
        image_plane_(left, right, bottom, top, distance, image_width,
                     image_height),
        tmo_(tmo) {
    w = -(gaze.normalize());
    if (left_handed) {
      u = w.cross(up.normalize()).normalize();
      v = u.cross(w).normalize();
    } else {
      u = up.normalize().cross(w).normalize();
      v = w.cross(u).normalize();
    }

    top_left_corner = e - w * image_plane_.distance + image_plane_.left * u +
                      image_plane_.top * v;
    s_u_constant =
        ((image_plane_.right - image_plane_.left) / image_plane_.width) * u;
    s_v_constant =
        ((image_plane_.top - image_plane_.bottom) / image_plane_.height) * v;
  }
  Camera(Camera&& rhs)
      : u(rhs.u),
        v(rhs.v),
        w(rhs.w),
        e(rhs.e),
        s_u_constant(rhs.s_u_constant),
        s_v_constant(rhs.s_v_constant),
        top_left_corner(rhs.top_left_corner),
        number_of_samples_(rhs.number_of_samples_),
        filename_(rhs.filename_),
        image_plane_(rhs.image_plane_),
        tmo_(rhs.tmo_) {
    rhs.tmo_ = nullptr;
  }

  ~Camera() {
    if (tmo_) {
      delete tmo_;
    }
  }
  Ray calculate_ray_at(float x, float y, float time = 0.0f) const {
    // For now, origin is top_left, right handed camera
    const Vector3 s = top_left_corner + x * s_u_constant - y * s_v_constant;
    Ray ray(e, (s - e).normalize());
    // ray.bg_u = x / image_plane_.width;
    // ray.bg_v = y / image_plane_.height;
    return ray;
  }

  const Image_plane& get_image_plane() const { return image_plane_; }
  const std::string& get_filename() const { return filename_; }
  int get_number_of_samples() const { return number_of_samples_; }
  const Tonemapping_operator* get_tmo() const { return tmo_; }

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
  Tonemapping_operator* tmo_;
};
#endif
