#ifndef TEXTURE_H
#define TEXTURE_H
//#define cimg_use_jpeg
#include "CImg.h"
#include "Vector3.h"
#include <string>
class Texture {
public:
  enum Interpolation_type {
    it_nearest,
    it_bilinear
  };
  enum Decal_mode {
    dm_replace_kd,
    dm_blend_kd,
    dm_replace_all
  };
  enum Appearance {
    a_repeat,
    a_clamp
  };
  Texture(const std::string& image_name, const std::string& interpolation_type, const std::string& decal_mode, const std::string& appearance)
  : texture_image(image_name.c_str()) {
    width_ = texture_image.width();
    height_ = texture_image.height();
    interpolation_type_ = to_interpolation_type(interpolation_type);
    decal_mode_ = to_decal_mode(decal_mode);
    appearance_ = to_appearance(appearance);
  }
  Vector3 get_color_at(float u, float v) const {
    if (appearance_ == a_clamp) {
      u = fmax(0.0f, fmin(1.0f, u));
      v = fmax(0.0f, fmin(1.0f, v));
    }
    else {
      u -= floor(u);
      v -= floor(v);
    }

    u *= width_;
    if (u >= width_) u--;
    v *= height_;
    if (v >= height_) v--;
    Vector3 color;
    if (interpolation_type_ == it_nearest) {
      color.x = texture_image((unsigned int) u, (unsigned int) v, 0);
      color.y = texture_image((unsigned int) u, (unsigned int) v, 1);
      color.z = texture_image((unsigned int) u, (unsigned int) v, 2);
    }
    else {
      const unsigned int p = u;
      const unsigned int q = v;
      const float dx = u - p;
      const float dy = v - q;
      color.x = texture_image(p, q, 0) * (1 - dx) * (1 - dy);
      color.x += texture_image(p+1, q, 0) * (dx) * (1 - dy);
      color.x += texture_image(p+1, q+1, 0) * (dx) * (dy);
      color.x += texture_image(p, q+1, 0) * (1 - dx) * (dy);

      color.y = texture_image(p, q, 1) * (1 - dx) * (1 - dy);
      color.y += texture_image(p + 1, q, 1) * (dx) * (1 - dy);
      color.y += texture_image(p + 1, q + 1, 1) * (dx) * (dy);
      color.y += texture_image(p, q + 1, 1) * (1 - dx) * (dy);

      color.z = texture_image(p, q, 2) * (1 - dx) * (1 - dy);
      color.z += texture_image(p + 1, q, 2) * (dx) * (1 - dy);
      color.z += texture_image(p + 1, q + 1, 2) * (dx) * (dy);
      color.z += texture_image(p, q + 1, 2) * (1 - dx) * (dy);
    }

    if (decal_mode_ == dm_blend_kd || decal_mode_ == dm_replace_kd) {
      color /= 255.;
    }
    return color;
  }

private:
  Interpolation_type interpolation_type_;
  Decal_mode decal_mode_;
  Appearance appearance_;
  cimg_library::CImg<unsigned char> texture_image;
  int width_;
  int height_;
  Interpolation_type to_interpolation_type(const std::string& str) {
    return str == "nearest" ? it_nearest : it_bilinear;
  }

  Decal_mode to_decal_mode(const std::string& str) {
    if (str == "replace_kd") return dm_replace_kd;
    return str == "blend_kd" ? dm_blend_kd : dm_replace_all;
  }

  Appearance to_appearance(const std::string& str) {
    return str == "repeat" ? a_repeat : a_clamp;
  }
};
#endif //TEXTURE_H
