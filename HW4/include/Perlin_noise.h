#ifndef PERLIN_NOISE_H_
#define PERLIN_NOISE_H_
#include "Vector3.h"
#include <vector>

class Perlin_noise {
public:
  enum Perlin_noise_appearance { pn_vein, pn_patch };
  Perlin_noise(const std::string& pn_appearance, float scaling_factor = 1.0f);
  float get_value_at(const Vector3& position) const;
private:
  float scaling_factor_;
  Perlin_noise_appearance pn_appearance_;
  std::vector<int> p_;
  static const std::vector<Vector3> g_;
  Vector3 edge_vector_at(int i, int j, int k) const;
  inline float weight_function(float x) const;
};
#endif
