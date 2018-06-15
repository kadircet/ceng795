#include "Light_sphere.h"
#include <random>
Vector3 Light_sphere::direction_and_distance(const Vector3& from_point,
                                             const Vector3& normal,
                                             float& distance,
                                             float& probability) const {
  thread_local static std::random_device rd;
  thread_local static std::mt19937 generator(rd());
  std::uniform_real_distribution<float> uniform_dist(0.0f, 1.0f);
  float epsilon_1 = uniform_dist(generator);
  float epsilon_2 = uniform_dist(generator);

  Vector3 point_in_sphere_space =
      transformation_.get_inverse_transformation_matrix().multiply(from_point);
  Vector3 w = center - point_in_sphere_space;
  float d = w.length();
  w = w.normalize();
  const Vector3 u = ((w.x != 0.0f || w.y != 0.0f) ? Vector3(-w.y, w.x, 0.0f)
                                                  : Vector3(0.0f, 1.0f, 0.0f))
                        .normalize();
  const Vector3 v = w.cross(u);
  float sin_theta_max = std::max(-1.0f, std::min(1.0f, radius / d));
  float theta_max = std::asin(sin_theta_max);
  float cos_theta_max = std::cos(theta_max);
  float phi = 2 * M_PI * epsilon_1;
  float theta = std::acos(std::max(
      -1.0f, std::min(1.0f, 1 - epsilon_2 + epsilon_2 * cos_theta_max)));

  Vector3 l_in_object_space =
      (w * std::cos(theta) + v * std::sin(theta) * std::cos(phi) +
       u * std::sin(theta) * std::sin(phi))
          .normalize();

  Vector3 l_in_world_space = transformation_.get_transformation_matrix()
                                 .multiply(l_in_object_space, true)
                                 .normalize();

  Ray ray_in_world_space(
      from_point + scene_->shadow_ray_epsilon * l_in_world_space,
      l_in_world_space, r_primary);
  Hit_data light_hit_data;
  light_hit_data.t = std::numeric_limits<float>::infinity();
  light_hit_data.shape = NULL;
  Sphere::intersect(ray_in_world_space, light_hit_data);
  // if (!light_hit_data.shape) {
  //  std::cout << "wtf" << std::endl;
  //}
  distance = light_hit_data.t;
  probability = 1 / (2 * M_PI * (1 - cos_theta_max));
  if (isnan(probability)) {
    std::cout << "wtf nan" << std::endl;
  }
  return l_in_world_space;
}

// Incoming radiance to the point from the light
Vector3 Light_sphere::incoming_radiance(const Vector3& from_point_to_light,
                                        float probability) const {
  return radiance_ / probability;
}
