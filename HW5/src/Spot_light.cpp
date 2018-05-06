#include "Spot_light.h"

Spot_light::Spot_light(const Vector3& position, const Vector3& intensity,
                       const Vector3& direction,
                       float coverage_angle_in_radians,
                       float falloff_angle_in_radians)
    : position_(position),
      intensity_(intensity),
      direction_(direction.normalize()) {
  cos_half_of_coverage_angle_ = std::cos(coverage_angle_in_radians / 2);
  cos_half_of_falloff_angle_ = std::cos(falloff_angle_in_radians / 2);
}

Vector3 Spot_light::direction_and_distance(const Vector3& from_point,
                                           float& distance) const {
  const Vector3 direction = position_ - from_point;
  distance = direction.length();
  return direction;
}

Vector3 Spot_light::incoming_radiance(
    const Vector3& from_point_to_light) const {
  Vector3 reverse_w_i = -(from_point_to_light.normalize());
  float cos_theta = reverse_w_i.dot(direction_);
  if (cos_theta > cos_half_of_falloff_angle_) {
    float x = from_point_to_light.x;
    float y = from_point_to_light.y;
    float z = from_point_to_light.z;
    return intensity_ / (x * x + y * y + z * z);
  } else if (cos_theta < cos_half_of_coverage_angle_) {
    return Vector3(0.0f);
  } else {
    float c = (cos_theta - cos_half_of_coverage_angle_) /
              (cos_half_of_falloff_angle_ - cos_half_of_coverage_angle_);
    float x = from_point_to_light.x;
    float y = from_point_to_light.y;
    float z = from_point_to_light.z;
    return intensity_ * pow(c, 4) / (x * x + y * y + z * z);
  }
  return Vector3();
}
