#ifndef BLINN_PHONG_BRDF_H_
#define BLINN_PHONG_BRDF_H_
#include "BRDF.h"
class Blinn_phong_BRDF : public BRDF {
 public:
  Blinn_phong_BRDF(float phong_exponent) : phong_exponent_(phong_exponent) {}
  Vector3 get_reflectance(const Vector3& normal, const Vector3& diffuse,
                          const Vector3& specular, const Vector3& w_i,
                          const Vector3& w_o) const override;

 private:
  float phong_exponent_;
};
#endif
