#include "Perlin_noise.h"
#include <chrono>
#include <cmath>
#include <random>

const std::vector<Vector3> Perlin_noise::g_ = {
    Vector3(1, 1, 0), Vector3(-1, 1, 0), Vector3(1, -1, 0), Vector3(-1, -1, 0),
    Vector3(1, 0, 1), Vector3(-1, 0, 1), Vector3(1, 0, -1), Vector3(-1, 0, -1),
    Vector3(0, 1, 1), Vector3(0, -1, 1), Vector3(0, 1, -1), Vector3(0, -1, -1),
    Vector3(1, 1, 0), Vector3(-1, 1, 0), Vector3(0, -1, 1), Vector3(0, -1, -1)};

Perlin_noise::Perlin_noise(const std::string& pn_appearance,
                           float scaling_factor)
    : scaling_factor_(scaling_factor) {
  pn_appearance_ = (pn_appearance == "vein") ? pn_vein : pn_patch;
  p_ = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  std::mt19937 random_engine;
  random_engine.seed(
      std::chrono::system_clock::now().time_since_epoch().count());
  std::shuffle(p_.begin(), p_.end(), random_engine);
}

float Perlin_noise::get_value_at(const Vector3& position) const {
  Vector3 scaled_position = position * scaling_factor_;
  float value = 0.0f;
  int x = floor(scaled_position.x);
  int y = floor(scaled_position.y);
  int z = floor(scaled_position.z);
  for (int i = x; i <= x + 1; i++) {
    for (int j = y; j <= y + 1; j++) {
      for (int k = z; k <= z + 1; k++) {
        const Vector3 v(scaled_position.x - i, scaled_position.y - j,
                        scaled_position.z - k);
        const Vector3 e = edge_vector_at(i, j, k);
        const float w_x = weight_function(v.x);
        const float w_y = weight_function(v.y);
        const float w_z = weight_function(v.z);
        value += e.dot(v) * w_x * w_y * w_z;
      }
    }
  }
  switch (pn_appearance_) {
    case pn_vein:
      return fabs(value);
    case pn_patch:
      return (value + 1.0f) / 2.0f;
  }
  return 0.0f;
}
inline float Perlin_noise::weight_function(float x) const {
  return -6 * pow(fabs(x), 5) + 15 * pow(fabs(x), 4) + -10 * pow(fabs(x), 3) +
         1;
}
Vector3 Perlin_noise::edge_vector_at(int i, int j, int k) const {
  i = i % 16;
  j = j % 16;
  k = k % 16;
  if (i < 0) i += 16;
  if (j < 0) j += 16;
  if (k < 0) k += 16;
  return g_[p_[(i + p_[(j + p_[k % 16]) % 16]) % 16]];
}
