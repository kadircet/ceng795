#include "Blinn_phong_BRDF.h"
#include <algorithm>
#include "Hit_data.h"
Vector3 Blinn_phong_BRDF::get_reflectance(const Hit_data& hit_data,
                                          const Vector3& diffuse,
                                          const Vector3& specular,
                                          const Vector3& w_i,
                                          const Vector3& w_o) const {
  const Vector3& normal = hit_data.normal;
  float cos_theta_i = normal.dot(w_i);
  if (cos_theta_i > 1.0f || cos_theta_i <= 0.0f) {
    return 0.0f;
  }
  Vector3 w_h = (w_i + w_o).normalize();
  float cos_alpha_h = w_h.dot(normal);
  Vector3 k_diffuse = diffuse;
  Vector3 k_specular =
      (specular * std::pow(cos_alpha_h, phong_exponent_)) / cos_theta_i;
  return k_specular + k_diffuse;
}
