#include "Light_sphere.h"
Vector3 Light_sphere::direction_and_distance(const Vector3& from_point,
                                             float& distance,
                                             float& probability) const {
  distance = 1.0f;
  probability = 1.0f;
  return 0.0f;
}

// Incoming radiance to the point from the light
Vector3 Light_sphere::incoming_radiance(const Vector3& from_point_to_light,
                                        float probability) const {
  return 0.0f;
}
