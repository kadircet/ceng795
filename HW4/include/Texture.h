#ifndef TEXTURE_H
#define TEXTURE_H
#include <string>
#include "Vector3.h"
#include "Perlin_noise.h"
class Texture {
 public:
  Texture(const std::string& image_name, const std::string& interpolation_type,
          const std::string& decal_mode, const std::string& appearance, const float normalizer, const float scaling_factor);
  Texture(Texture&& rhs);
  ~Texture();
  Vector3 get_color_at(float u, float v) const;
  Vector3 blend_color(const Vector3& texture_color, const Vector3& diffuse_color) const;
  enum Interpolation_type { it_nearest, it_bilinear };
  enum Decal_mode { dm_replace_kd, dm_blend_kd, dm_replace_all };
  enum Appearance { a_repeat, a_clamp };

  inline Decal_mode get_decal_mode() const { return decal_mode_; }
  inline Interpolation_type get_interpolation_type() const {
    return interpolation_type_;
  }
  inline Appearance get_appearance_() const { return appearance_; }
  inline const Perlin_noise* get_perlin_noise() const { return perlin_noise_;}
  inline bool is_perlin_noise() const {
    return perlin_noise_ != nullptr;
  }
 private:
  Interpolation_type interpolation_type_;
  Decal_mode decal_mode_;
  Appearance appearance_;
  unsigned char* texture_image_;
  int width_;
  int height_;
  float normalizer_;
  Perlin_noise* perlin_noise_;
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
#endif  // TEXTURE_H
