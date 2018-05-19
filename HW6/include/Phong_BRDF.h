#ifndef PHONG_BRDF_H_
#define PHONG_BRDF_H_
#include "BRDF.h"
class Phong_BRDF : public BRDF {
 public:
  Phong_BRDF(float phong_exponent, bool normalized)
      : phong_exponent_(phong_exponent), normalized_(normalized) {}
  Vector3 get_reflectance(const Hit_data& hit_data, const Vector3& diffuse,
                          const Vector3& specular, const Vector3& w_i,
                          const Vector3& w_o) const override;

 private:
  float phong_exponent_;
  bool normalized_;
};
#endif
