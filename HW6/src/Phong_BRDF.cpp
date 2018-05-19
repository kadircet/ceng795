#include "Phong_BRDF.h"
#include <algorithm>
#include "Hit_data.h"
Vector3 Phong_BRDF::get_reflectance(const Hit_data& hit_data,
                                    const Vector3& diffuse,
                                    const Vector3& specular, const Vector3& w_i,
                                    const Vector3& w_o) const {
  const Vector3& normal = hit_data.normal;
  float theta_i = normal.dot(w_i);
  Vector3 reflected_wi = 2 * theta_i * normal - w_i;
  float alpha = std::max(0.0f, w_o.dot(reflected_wi));
  Vector3 k_diffuse = diffuse / M_PI;
  Vector3 k_specular = (specular * std::pow(alpha, phong_exponent_)) / theta_i;
  if (normalized_) {
    k_specular = ((k_specular) * (phong_exponent_ + 1.0f)) / (2.0f * M_PI);
  }
  return (k_specular + k_diffuse) * theta_i;
}
