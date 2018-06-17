#include "Mesh_triangle.h"
#include <algorithm>
#include "Intersection.h"
#include "Scene.h"
Mesh_triangle::Mesh_triangle(const Scene* scene, int vertex_index_0,
                             int vertex_index_1, int vertex_index_2,
                             int vertex_offset, int texture_offset,
                             int material_id, int texture_id,
                             Triangle_shading_mode tsm)
    : vertex_index_0(vertex_index_0),
      vertex_index_1(vertex_index_1),
      vertex_index_2(vertex_index_2),
      vertex_offset(vertex_offset),
      texture_offset(texture_offset),
      material_id(material_id),
      texture_id(texture_id),
      triangle_shading_mode(tsm),
      scene_(scene) {
  const Vector3& v_0 = scene_->get_vertex_at(vertex_index_0 + vertex_offset)
                           .get_vertex_position();
  const Vector3& v_1 = scene_->get_vertex_at(vertex_index_1 + vertex_offset)
                           .get_vertex_position();
  const Vector3& v_2 = scene_->get_vertex_at(vertex_index_2 + vertex_offset)
                           .get_vertex_position();
  normal = (v_1 - v_0).cross(v_2 - v_0).normalize();
  Vector3 min_c = v_0;
  Vector3 max_c = v_0;

  min_c.x = std::min(min_c.x, v_1.x);
  min_c.y = std::min(min_c.y, v_1.y);
  min_c.z = std::min(min_c.z, v_1.z);

  max_c.x = std::max(max_c.x, v_1.x);
  max_c.y = std::max(max_c.y, v_1.y);
  max_c.z = std::max(max_c.z, v_1.z);

  min_c.x = std::min(min_c.x, v_2.x);
  min_c.y = std::min(min_c.y, v_2.y);
  min_c.z = std::min(min_c.z, v_2.z);

  max_c.x = std::max(max_c.x, v_2.x);
  max_c.y = std::max(max_c.y, v_2.y);
  max_c.z = std::max(max_c.z, v_2.z);
  bounding_box_ = Bounding_box(min_c, max_c);
}

float Mesh_triangle::get_surface_area() const {
  const Vector3& v_0 = scene_->get_vertex_at(vertex_index_0 + vertex_offset)
                           .get_vertex_position();
  const Vector3& v_1 = scene_->get_vertex_at(vertex_index_1 + vertex_offset)
                           .get_vertex_position();
  const Vector3& v_2 = scene_->get_vertex_at(vertex_index_2 + vertex_offset)
                           .get_vertex_position();
  return (v_1 - v_0).cross(v_2 - v_0).length() / 2;
}
bool Mesh_triangle::intersect(const Ray& ray, Intersection& intersection,
                              bool culling) const {
  const Vertex& vertex_0 =
      scene_->get_vertex_at(vertex_index_0 + vertex_offset);
  const Vertex& vertex_1 =
      scene_->get_vertex_at(vertex_index_1 + vertex_offset);
  const Vertex& vertex_2 =
      scene_->get_vertex_at(vertex_index_2 + vertex_offset);
  const Vector3& p_0 = vertex_0.get_vertex_position();
  const Vector3& p_1 = vertex_1.get_vertex_position();
  const Vector3& p_2 = vertex_2.get_vertex_position();
  const Vector3 a_col1 = p_0 - p_1;
  const Vector3 a_col2 = p_0 - p_2;
  const Vector3& a_col3 = ray.d;

  if (culling && a_col3.dot(normal) > 0.0f) {
    return false;
  }

  const float det_a = determinant(a_col1, a_col2, a_col3);
  if (det_a == 0.0f) {
    return false;
  }
  const Vector3 b = (p_0 - ray.o) / det_a;
  const float beta = determinant(b, a_col2, a_col3);
  if (beta < -intersection_test_epsilon) return false;
  const float gamma = determinant(a_col1, b, a_col3);
  if (gamma < -intersection_test_epsilon ||
      beta + gamma > 1.0f + intersection_test_epsilon) {
    return false;
  }
  const float t = determinant(a_col1, a_col2, b);
  if (t > -intersection_test_epsilon) {
    intersection.t = t;
    Vector3 local_intersection_point = ray.point_at(t);
    intersection.shape = this;
    // float u = -1;
    // float v = -1;
    // float perlin_value = -1;
    Vector3 normal;
    switch (triangle_shading_mode) {
      case tsm_smooth:
        normal = ((1 - beta - gamma) * vertex_0.get_vertex_normal() +
                  beta * vertex_1.get_vertex_normal() +
                  gamma * vertex_2.get_vertex_normal())
                     .normalize();
        break;
      case tsm_flat:
        normal = this->normal;
        break;
    }
    intersection.normal = normal;
    return true;
  }
  return false;
}
