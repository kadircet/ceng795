#include "Torrance_sparrow_BRDF.h"
#include <algorithm>
#include "Hit_data.h"
Vector3 Torrance_sparrow_BRDF::get_reflectance(const Hit_data& hit_data,
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
  float cos_alpha = w_h.dot(normal);
  float d_alpha =
      ((exponent_ + 2) * std::pow(cos_alpha, exponent_)) / (2.0f * M_PI);
  float cos_theta_o = w_o.dot(normal);
  float cos_beta = w_o.dot(w_h);
  float g_wi_wo = std::min(2 * cos_alpha * cos_theta_o / cos_beta,
                           2 * cos_alpha * cos_theta_i / cos_beta);
  g_wi_wo = g_wi_wo = std::min(1.0f, g_wi_wo);
  float r_0 = std::pow(refractive_index_ - 1.0f, 2) /
              std::pow(refractive_index_ + 1.0f, 2);
  float f_beta = r_0 + (1 - r_0) * std::pow(1 - cos_beta, 5);
  Vector3 k_diffuse = diffuse / M_PI;
  Vector3 k_specular =
      specular * d_alpha * f_beta * g_wi_wo / (4 * cos_theta_i * cos_theta_o);
  return k_diffuse + k_specular;
}
