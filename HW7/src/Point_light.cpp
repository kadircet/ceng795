#include "Point_light.h"

Point_light::Point_light(const Vector3& position, const Vector3& intensity)
    : position_(position), intensity_(intensity) {}

Vector3 Point_light::direction_and_distance(
    const Vector3& from_point, float& distance,
    float& probability) const {
  const Vector3 direction = position_ - from_point;
  distance = direction.length();
  return direction;
}

Vector3 Point_light::incoming_radiance(const Vector3& from_point_to_light,
                                       float probability) const {
  float x = from_point_to_light.x;
  float y = from_point_to_light.y;
  float z = from_point_to_light.z;
  return intensity_ / (x * x + y * y + z * z);
}
