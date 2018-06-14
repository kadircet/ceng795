#include "Spherical_directional_light.h"
#include <random>
#include "tinyexr.h"
Spherical_directional_light::Spherical_directional_light(
    const std::string& envmap_name) {
  const char* err = "error";
  int ret = LoadEXR(&env_map_, &width_, &height_, envmap_name.c_str(), &err);
  if (ret != 0) {
    std::cerr << "Cannot load EXR file." << std::endl;
    exit(-1);
  };
  // TODO free env_map_
}

Vector3 Spherical_directional_light::direction_and_distance(
    const Vector3& from_point, const Vector3& normal, float& distance,
    float& probability) const {
  thread_local static std::random_device rd;
  thread_local static std::mt19937 generator(rd());
  std::uniform_real_distribution<float> uniform_dist(0.0f, 1.0f);
  float epsilon_1 = uniform_dist(generator);
  float epsilon_2 = uniform_dist(generator);

  Vector3 w = normal;
  const Vector3 u = ((w.x != 0.0f || w.y != 0.0f) ? Vector3(-w.y, w.x, 0.0f)
                                                  : Vector3(0.0f, 1.0f, 0.0f))
                        .normalize();
  const Vector3 v = w.cross(u);
  float phi = 2 * M_PI * epsilon_1;
  float theta = std::acos(epsilon_2);
  distance = std::numeric_limits<float>::max();
  probability = 1.0f / (2.0f * M_PI);
  return (w * std::cos(theta) + v * std::sin(theta) * std::cos(phi) +
          u * std::sin(theta) * std::sin(phi))
      .normalize();
  ;
}

Vector3 Spherical_directional_light::incoming_radiance(
    const Vector3& from_point_to_light, float probability) const {
  float cos_theta = from_point_to_light.y;
  float theta = std::acos(clamp(-1.0f, 1.0f, cos_theta));
  float phi = std::atan2(from_point_to_light.z, from_point_to_light.x);
  float u = (M_PI - phi) / (2.0f * M_PI);
  float v = theta / M_PI;
  u = std::max(0.0f, std::min(1.0f, u));
  v = std::max(0.0f, std::min(1.0f, v));
  u *= width_;
  if (u >= width_) u--;
  v *= height_;
  if (v >= height_) v--;
  Vector3 radiance;
  unsigned int coord = 4 * width_ * (unsigned int)v + 4 * (unsigned int)u;
  radiance.x = env_map_[coord];
  radiance.y = env_map_[coord + 1];
  radiance.z = env_map_[coord + 2];
  return radiance / probability;
}
