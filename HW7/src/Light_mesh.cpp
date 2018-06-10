#include "Light_mesh.h"
#include <random>
Vector3 Light_mesh::direction_and_distance(const Vector3& from_point,
                                           float& distance,
                                           float& probability) const {
  thread_local static std::random_device rd;
  thread_local static std::mt19937 generator(rd());
  std::uniform_real_distribution<float> uniform_dist(0.0f, 1.0f);
  float epsilon_0 = uniform_dist(generator);
  float epsilon_1 = std::sqrt(uniform_dist(generator));
  float epsilon_2 = uniform_dist(generator);
  auto x = std::lower_bound(triangles_.begin(), triangles_.end(),
                            std::pair<const float, Shape*>(epsilon_0, nullptr));
  Mesh_triangle* triangle = (Mesh_triangle*)x->second;
  const Vertex& vertex_0 =
      scene_->get_vertex_at(triangle->index_0 + triangle->offset);
  const Vertex& vertex_1 =
      scene_->get_vertex_at(triangle->index_0 + triangle->offset);
  const Vertex& vertex_2 =
      scene_->get_vertex_at(triangle->index_0 + triangle->offset);
  const Vector3& p_0 = vertex_0.get_vertex_position();
  const Vector3& p_1 = vertex_1.get_vertex_position();
  const Vector3& p_2 = vertex_2.get_vertex_position();

  Vector3 q = (1.0f - epsilon_2) * p_1 + epsilon_2 * p_2;
  Vector3 p = (1.0f - epsilon_1) * p_0 + epsilon_1 * q;

  const Vector3 direction = p - from_point;
  distance = direction.length();
  float cos_theta_i =
      std::max(0.001f, direction.normalize().dot(triangle->normal));
  probability = distance * distance / (total_area_ * cos_theta_i);
  return direction;
  ;
}

// Incoming radiance to the point from the light
Vector3 Light_mesh::incoming_radiance(const Vector3& from_point_to_light,
                                      float probability) const {
  return radiance_ * 255.0f / probability;
}
