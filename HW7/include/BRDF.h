#ifndef BRDF_H_
#define BRDF_H_
#include "Vector3.h"
class Hit_data;
class BRDF {
 public:
  virtual Vector3 get_reflectance(const Hit_data& hit_data,
                                  const Vector3& diffuse,
                                  const Vector3& specular, const Vector3& w_i,
                                  const Vector3& w_o) const = 0;
};
#endif
