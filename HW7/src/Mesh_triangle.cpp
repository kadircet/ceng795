#include "Mesh_triangle.h"
#include <algorithm>
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
  const Vector3& p_0 = vertex_0.get_vertex_position();
  const Vector3& p_1 = vertex_1.get_vertex_position();
  const Vector3& p_2 = vertex_2.get_vertex_position();
  const Vector3 a_col1 = p_0 - p_1;
  const Vector3 a_col2 = p_0 - p_2;
  const Vector3& a_col3 = ray.d;
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
    hit_data.t = t;
    Vector3 local_intersection_point = ray.point_at(t);
    hit_data.shape = this;
    float u = -1;
    float v = -1;
    float perlin_value = -1;
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
    if (texture_id != -1) {
      // TODO: Maybe do this calculations outside to give different textures to
      // different mesh instances like materials
      const Texture& texture = scene_->textures[texture_id];
      if (texture.is_perlin_noise()) {
        perlin_value =
            texture.get_perlin_noise()->get_value_at(local_intersection_point);
        if (texture.is_bump()) {
          normal = texture.bump_normal(normal, local_intersection_point);
        }
      } else {
        Vector3 uva = scene_->texture_coord_data[index_0 + texture_offset];
        Vector3 uvb = scene_->texture_coord_data[index_1 + texture_offset];
        Vector3 uvc = scene_->texture_coord_data[index_2 + texture_offset];
        u = uva.x + beta * (uvb.x - uva.x) + gamma * (uvc.x - uva.x);
        v = uva.y + beta * (uvb.y - uva.y) + gamma * (uvc.y - uva.y);
        if (texture.is_bump()) {
          // calculate gradients
          float ub_ua = uvb.x - uva.x;
          float uc_ua = uvc.x - uva.x;
          float vb_va = uvb.y - uva.y;
          float vc_va = uvc.y - uva.y;
          Vector3 pb_pa = p_1 - p_0;
          Vector3 pc_pa = p_2 - p_0;
          float inverse_constant = (ub_ua * vc_va) - (vb_va * uc_ua);
          if (inverse_constant == 0.0f) {
            std::cerr << "Inverse constant == 0 at triangle bump mapping"
                      << std::endl;
            inverse_constant = 0.000001f;
          }
          inverse_constant = 1.0f / inverse_constant;
          Vector3 dp_du = inverse_constant * (vc_va * pb_pa - vb_va * pc_pa);
          Vector3 dp_dv = inverse_constant * (ub_ua * pc_pa - uc_ua * pb_pa);
          // std::cerr<<normal<<std::endl<<std::endl;
          // std::cerr<<dp_dv.cross(dp_du).normalize()<<std::endl;
          normal = texture.bump_normal(dp_du, dp_dv, normal, u, v);
        }
      }
    }
    hit_data.u = u;
    hit_data.v = v;
    hit_data.perlin_value = perlin_value;
    hit_data.normal = normal;
    hit_data.is_light_object = false;
    return true;
  }
  return false;
}
