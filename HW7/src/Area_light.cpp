#include "Area_light.h"
#include <random>

Area_light::Area_light(const Vector3& position, const Vector3& intensity,
                       const Vector3& edge_vector_1,
                       const Vector3& edge_vector_2)
    : position_(position),
      intensity_(intensity),
      edge_vector_1_(edge_vector_1),
      edge_vector_2_(edge_vector_2) {
  normal_ = edge_vector_1_.cross(edge_vector_2_).normalize();
}

Vector3 Area_light::direction_and_distance(const Vector3& from_point,
                                           float& distance,
                                           float& probability) const {
  thread_local static std::random_device rd;
  thread_local static std::mt19937 generator(rd());
  std::uniform_real_distribution<float> area_light_distribution(0.0f, 1.0f);
  float epsilon_1 = area_light_distribution(generator);
  float epsilon_2 = area_light_distribution(generator);
  const Vector3 position =
      position_ + edge_vector_1_ * epsilon_1 + edge_vector_2_ * epsilon_2;
  const Vector3 direction = position - from_point;
  distance = direction.length();
  return direction;
}

Vector3 Area_light::incoming_radiance(const Vector3& from_point_to_light,
                                      float probability) const {
  Vector3 reverse_w_i = -(from_point_to_light.normalize());
  float x = from_point_to_light.x;
  float y = from_point_to_light.y;
  float z = from_point_to_light.z;
  return intensity_ * (reverse_w_i.dot(normal_)) / (x * x + y * y + z * z);
}
