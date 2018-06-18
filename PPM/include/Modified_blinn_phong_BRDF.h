#ifndef MODIFIED_BLINN_PHONG_BRDF_H_
#define MODIFIED_BLINN_PHONG_BRDF_H_
#include "BRDF.h"
class Modified_blinn_phong_BRDF : public BRDF {
 public:
  Modified_blinn_phong_BRDF(float phong_exponent, bool normalized)
      : phong_exponent_(phong_exponent), normalized_(normalized) {}
  Vector3 get_reflectance(const Vector3& normal, const Vector3& diffuse,
                          const Vector3& specular, const Vector3& w_i,
                          const Vector3& w_o) const override;

 private:
  float phong_exponent_;
  bool normalized_;
};
#endif
