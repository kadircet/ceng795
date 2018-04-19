#include "Texture.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

Texture::Texture(const std::string& image_name,
                 const std::string& interpolation_type,
                 const std::string& decal_mode, const std::string& appearance) {
  int n;
  texture_image_ = stbi_load(image_name.c_str(), &width_, &height_, &n, 3);
  interpolation_type_ = to_interpolation_type(interpolation_type);
  decal_mode_ = to_decal_mode(decal_mode);
  appearance_ = to_appearance(appearance);
}
Texture::Texture(Texture&& rhs)
    : texture_image_(rhs.texture_image_),
      width_(rhs.width_),
      height_(rhs.height_),
      interpolation_type_(rhs.interpolation_type_),
      decal_mode_(rhs.decal_mode_),
      appearance_(rhs.appearance_) {
  rhs.texture_image_ = nullptr;
}

Texture::~Texture() {
  if (texture_image_) stbi_image_free(texture_image_);
}

Vector3 Texture::get_color_at(float u, float v) const {
  if (appearance_ == a_clamp) {
    u = fmax(0.0f, fmin(1.0f, u));
    v = fmax(0.0f, fmin(1.0f, v));
  } else {
    u = fmax(0.0f, fmin(1.0f, u - (unsigned int)u));
    v = fmax(0.0f, fmin(1.0f, v - (unsigned int)v));
  }

  u *= width_;
  // if (u >= width_) u--;
  v *= height_;
  // if (v >= height_) v--;
  Vector3 color;
  if (interpolation_type_ == it_nearest) {
    unsigned int coord = 3 * width_ * (unsigned int)v + 3 * (unsigned int)u;
    color.x = texture_image_[coord];
    color.y = texture_image_[coord + 1];
    color.z = texture_image_[coord + 2];
  } else {
    const unsigned int p = u;
    const unsigned int q = v;
    const float dx = u - p;
    const float dy = v - q;
    const int c_p_q = 3 * width_ * q + 3 * p;
    const int c_pn_q = c_p_q + 3;
    const int c_p_qn = c_p_q + 3 * width_;
    const int c_pn_qn = c_p_qn + 3;
    color.x = texture_image_[c_p_q] * (1 - dx) * (1 - dy);
    color.x += texture_image_[c_pn_q] * (dx) * (1 - dy);
    color.x += texture_image_[c_pn_qn] * (dx) * (dy);
    color.x += texture_image_[c_p_qn] * (1 - dx) * (dy);

    color.y = texture_image_[c_p_q + 1] * (1 - dx) * (1 - dy);
    color.y += texture_image_[c_pn_q + 1] * (dx) * (1 - dy);
    color.y += texture_image_[c_pn_qn + 1] * (dx) * (dy);
    color.y += texture_image_[c_p_qn + 1] * (1 - dx) * (dy);

    color.z = texture_image_[c_p_q + 2] * (1 - dx) * (1 - dy);
    color.z += texture_image_[c_pn_q + 2] * (dx) * (1 - dy);
    color.z += texture_image_[c_pn_qn + 2] * (dx) * (dy);
    color.z += texture_image_[c_p_qn + 2] * (1 - dx) * (dy);
  }

  if (decal_mode_ == dm_blend_kd || decal_mode_ == dm_replace_kd) {
    color /= 255.0f;
  }
  return color;
}
