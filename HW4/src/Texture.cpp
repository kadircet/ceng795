#include "Texture.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

Texture::Texture(const std::string& image_name,
                 const std::string& interpolation_type,
                 const std::string& decal_mode, const std::string& appearance,
                 const float normalizer, const float scaling_factor,
                 const bool is_bump, const float bumpmap_multiplier) {
  std::cout << decal_mode << std::endl;
  int n;
  if (image_name != std::string("perlin")) {
    texture_image_ = stbi_load(image_name.c_str(), &width_, &height_, &n, 3);
    perlin_noise_ = nullptr;
  } else {
    texture_image_ = nullptr;
    perlin_noise_ = new Perlin_noise(appearance, scaling_factor);
  }
  appearance_ = to_appearance(appearance);
  interpolation_type_ = to_interpolation_type(interpolation_type);
  decal_mode_ = to_decal_mode(decal_mode);
  normalizer_ = normalizer;
  is_bump_ = is_bump;
  bumpmap_multiplier_ = bumpmap_multiplier;
}
Texture::Texture(Texture&& rhs)
    : interpolation_type_(rhs.interpolation_type_),
      decal_mode_(rhs.decal_mode_),
      appearance_(rhs.appearance_),
      texture_image_(rhs.texture_image_),
      width_(rhs.width_),
      height_(rhs.height_),
      normalizer_(rhs.normalizer_),
      perlin_noise_(rhs.perlin_noise_),
      is_bump_(rhs.is_bump_),
      bumpmap_multiplier_(rhs.bumpmap_multiplier_) {
  rhs.texture_image_ = nullptr;
  rhs.perlin_noise_ = nullptr;
}

Texture::~Texture() {
  if (texture_image_) stbi_image_free(texture_image_);
  if (perlin_noise_) free(perlin_noise_);
}
Vector3 Texture::blend_color(const Vector3& texture_color,
                             const Vector3& diffuse_color) const {
  switch (decal_mode_) {
    case dm_replace_kd:
      return texture_color;
    case dm_blend_kd:
      return (texture_color + diffuse_color) / 2.0f;
    case dm_replace_all:
      return Vector3(0.0f, 0.0f, 0.0f);
    default:
      return Vector3(0.0f, 0.0f, 0.0f);
  }
}
Vector3 Texture::get_color_at(float u, float v) const {
  if (appearance_ == a_clamp) {
    u = fmax(0.0f, fmin(1.0f, u));
    v = fmax(0.0f, fmin(1.0f, v));
  } else {
    u = fmax(0.0f, fmin(1.0f, u - floor(u)));
    v = fmax(0.0f, fmin(1.0f, v - floor(v)));
  }
  u *= width_;
  if (u >= width_) u--;
  v *= height_;
  if (v >= height_) v--;
  Vector3 color;
  if (interpolation_type_ == it_nearest) {
    unsigned int coord = 3 * width_ * (unsigned int)v + 3 * (unsigned int)u;
    color.x = texture_image_[coord];
    color.y = texture_image_[coord + 1];
    color.z = texture_image_[coord + 2];
  } else {
    if (appearance_ == a_repeat) {
      const unsigned int p = ((unsigned int)u) % width_;
      const unsigned int pn = (p + 1) % width_;
      const unsigned int q = ((unsigned int)v) % height_;
      const unsigned int qn = (q + 1) % height_;
      const float dx = u - p;
      const float dy = v - q;
      const int c_p_q = 3 * (width_ * q + p);
      const int c_pn_q = 3 * (width_ * q + pn);
      const int c_p_qn = 3 * (width_ * qn + p);
      const int c_pn_qn = 3 * (width_ * qn + pn);
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
    } else {
      const unsigned int p = u;
      const unsigned int pn = p + 1;
      const unsigned int q = v;
      const unsigned int qn = v + 1;
      const float dx = u - p;
      const float dy = v - q;
      const int c_p_q = 3 * (width_ * q + p);
      const int c_pn_q = 3 * (width_ * q + pn);
      const int c_p_qn = 3 * (width_ * qn + p);
      const int c_pn_qn = 3 * (width_ * qn + pn);
      color.x = texture_image_[c_p_q] * (1 - dx) * (1 - dy);
      color.y = texture_image_[c_p_q + 1] * (1 - dx) * (1 - dy);
      color.z = texture_image_[c_p_q + 2] * (1 - dx) * (1 - dy);
      if (pn < width_) {
        color.x += texture_image_[c_pn_q] * (dx) * (1 - dy);
        color.y += texture_image_[c_pn_q + 1] * (dx) * (1 - dy);
        color.z += texture_image_[c_pn_q + 2] * (dx) * (1 - dy);
      }
      if (qn < height_) {
        color.x += texture_image_[c_pn_qn] * (dx) * (dy);
        color.y += texture_image_[c_pn_qn + 1] * (dx) * (dy);
        color.z += texture_image_[c_pn_qn + 2] * (dx) * (dy);
      }
      if (pn < width_ && qn < height_) {
        color.x += texture_image_[c_p_qn] * (1 - dx) * (dy);
        color.y += texture_image_[c_p_qn + 1] * (1 - dx) * (dy);
        color.z += texture_image_[c_p_qn + 2] * (1 - dx) * (dy);
      }
    }
  }

  return color;
}

Vector3 Texture::bump_normal(const Vector3 &normal, const Vector3 &position) const {
  const float perlin_value = perlin_noise_->get_value_at(position);
  const float epsilon = 0.001;
  const float dd_dx = (perlin_noise_->get_value_at(position+Vector3(epsilon,0,0)) - perlin_value) / epsilon;
  const float dd_dy = (perlin_noise_->get_value_at(position+Vector3(0,epsilon,0)) - perlin_value) / epsilon;
  const float dd_dz = (perlin_noise_->get_value_at(position+Vector3(0,0,epsilon)) - perlin_value) / epsilon;
  Vector3 g = bumpmap_multiplier_*Vector3(dd_dx,dd_dy,dd_dz);
  Vector3 g_ii = normal*(g.dot(normal));
  Vector3 surface_g = g-g_ii;
  return normal-surface_g;
}


Vector3 Texture::bump_normal(const Vector3 &dp_du, const Vector3 &dp_dv, const Vector3 &normal, const float u, const float v) const {
  const Vector3 color_u_v = get_color_at(u, v);
  const Vector3 color_un_v = get_color_at(u+1.0f/width_, v);
  const Vector3 color_u_vn = get_color_at(u, v+1.0f/height_);
  const float value_u_v = (color_u_v.x+color_u_v.y+color_u_v.z) / 3.0f;
  const float value_un_v = (color_un_v.x+color_un_v.y+color_un_v.z) / 3.0f;
  const float value_u_vn = (color_u_vn.x+color_u_vn.y+color_u_vn.z) / 3.0f;
  
  const float dd_du = (value_un_v-value_u_v)*bumpmap_multiplier_;
  const float dd_dv = (value_u_vn-value_u_v)*bumpmap_multiplier_;
  const Vector3 dq_du = dp_du + dd_du*normal;
  const Vector3 dq_dv = dp_dv + dd_dv*normal;
  return dq_dv.cross(dq_du).normalize();
}

