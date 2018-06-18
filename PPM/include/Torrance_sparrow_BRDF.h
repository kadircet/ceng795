#ifndef TORRANCE_SPARROW_BRDF_H_
#define TORRANCE_SPARROW_BRDF_H_
#include "BRDF.h"
class Torrance_sparrow_BRDF : public BRDF {
 public:
  Torrance_sparrow_BRDF(float exponent, float refractive_index)
      : exponent_(exponent), refractive_index_(refractive_index) {}
  Vector3 get_reflectance(const Vector3& normal, const Vector3& diffuse,
                          const Vector3& specular, const Vector3& w_i,
                          const Vector3& w_o) const override;

 private:
  float refractive_index_;
  float exponent_;
};
#endif
