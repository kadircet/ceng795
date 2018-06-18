#ifndef MODIFIED_PHONG_BRDF_H_
#define MODIFIED_PHONG_BRDF_H_
#include "BRDF.h"
class Modified_phong_BRDF : public BRDF {
 public:
  Modified_phong_BRDF(float phong_exponent, bool normalized)
      : phong_exponent_(phong_exponent), normalized_(normalized) {}
  Vector3 get_reflectance(const Vector3& normal, const Vector3& diffuse,
                          const Vector3& specular, const Vector3& w_i,
                          const Vector3& w_o) const override;

 private:
  float phong_exponent_;
  bool normalized_;
};
#endif
