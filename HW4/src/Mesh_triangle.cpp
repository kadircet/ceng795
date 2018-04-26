#include "Mesh_triangle.h"
#include "Hit_data.h"
#include "Scene.h"
Mesh_triangle::Mesh_triangle(const Scene* scene, int index_0, int index_1,
                             int index_2, int offset, int texture_offset,
                             int material_id, int texture_id,
                             Triangle_shading_mode tsm)
    : index_0(index_0),
      index_1(index_1),
      index_2(index_2),
      offset(offset),
      texture_offset(texture_offset),
      material_id(material_id),
      texture_id(texture_id),
      triangle_shading_mode(tsm),
      scene_(scene) {
  const Vector3& v_0 =
      scene_->get_vertex_at(index_0 + offset).get_vertex_position();
  const Vector3& v_1 =
      scene_->get_vertex_at(index_1 + offset).get_vertex_position();
  const Vector3& v_2 =
      scene_->get_vertex_at(index_2 + offset).get_vertex_position();
  normal = (v_1 - v_0).cross(v_2 - v_0).normalize();
  Vector3 min_c = v_0;
  Vector3 max_c = v_0;

  min_c.x = fmin(min_c.x, v_1.x);
  min_c.y = fmin(min_c.y, v_1.y);
  min_c.z = fmin(min_c.z, v_1.z);

  max_c.x = fmax(max_c.x, v_1.x);
  max_c.y = fmax(max_c.y, v_1.y);
  max_c.z = fmax(max_c.z, v_1.z);

  min_c.x = fmin(min_c.x, v_2.x);
  min_c.y = fmin(min_c.y, v_2.y);
  min_c.z = fmin(min_c.z, v_2.z);

  max_c.x = fmax(max_c.x, v_2.x);
  max_c.y = fmax(max_c.y, v_2.y);
  max_c.z = fmax(max_c.z, v_2.z);
  bounding_box_ = Bounding_box(min_c, max_c);
}

float Mesh_triangle::get_surface_area() const {
  const Vector3& v_0 =
      scene_->get_vertex_at(index_0 + offset).get_vertex_position();
  const Vector3& v_1 =
      scene_->get_vertex_at(index_1 + offset).get_vertex_position();
  const Vector3& v_2 =
      scene_->get_vertex_at(index_2 + offset).get_vertex_position();
  return (v_1 - v_0).cross(v_2 - v_0).length() / 2;
}
bool Mesh_triangle::intersect(const Ray& ray, Hit_data& hit_data) const {
  const Vertex& vertex_0 = scene_->get_vertex_at(index_0 + offset);
  const Vertex& vertex_1 = scene_->get_vertex_at(index_1 + offset);
  const Vertex& vertex_2 = scene_->get_vertex_at(index_2 + offset);
  const Vector3& v_0 = vertex_0.get_vertex_position();
  const Vector3& v_1 = vertex_1.get_vertex_position();
  const Vector3& v_2 = vertex_2.get_vertex_position();
  const Vector3 a_col1 = v_0 - v_1;
  const Vector3 a_col2 = v_0 - v_2;
  const Vector3& a_col3 = ray.d;
  const float det_a = determinant(a_col1, a_col2, a_col3);
  if (det_a == 0.0f) {
    return false;
  }
  const Vector3 b = (v_0 - ray.o) / det_a;
  const float beta = determinant(b, a_col2, a_col3);
  if (beta < -intersection_test_epsilon) return false;
  const float gamma = determinant(a_col1, b, a_col3);
  if (gamma < -intersection_test_epsilon ||
      beta + gamma > 1.0f + intersection_test_epsilon) {
    return false;
  }
  const float t = determinant(a_col1, a_col2, b);
  if (t > -intersection_test_epsilon) {
    hit_data.t = t;
    hit_data.shape = this;
    if (texture_id != -1 && !scene_->textures[texture_id].is_perlin_noise()) {
      Vector3 ua = scene_->texture_coord_data[index_0 + texture_offset];
      Vector3 ub = scene_->texture_coord_data[index_1 + texture_offset];
      Vector3 uc = scene_->texture_coord_data[index_2 + texture_offset];
      hit_data.u = ua.x + beta * (ub.x - ua.x) + gamma * (uc.x - ua.x);
      hit_data.v = ua.y + beta * (ub.y - ua.y) + gamma * (uc.y - ua.y);
    } else {
      hit_data.u = -1;
      hit_data.v = -1;
    }
    switch (triangle_shading_mode) {
      case tsm_smooth:
        hit_data.normal = ((1 - beta - gamma) * vertex_0.get_vertex_normal() +
                           beta * vertex_1.get_vertex_normal() +
                           gamma * vertex_2.get_vertex_normal())
                              .normalize();
        break;
      case tsm_flat:
        hit_data.normal = this->normal;
        break;
    }
    return true;
  }
  return false;
}
