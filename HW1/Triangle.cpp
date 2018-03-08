#include "Triangle.h"
#include "Hit_data.h"
#include "Scene.h"
Triangle::Triangle(const Scene* scene, int index_0, int index_1, int index_2,
                   int material_id)
    : index_0(index_0),
      index_1(index_1),
      index_2(index_2),
      material_id(material_id),
      scene_(scene) {
  const Vector3& v_0 = scene_->get_vertex_at(index_0);
  const Vector3& v_1 = scene_->get_vertex_at(index_1);
  const Vector3& v_2 = scene_->get_vertex_at(index_2);
  normal = (v_1 - v_0).cross(v_2 - v_0).normalize();
}
Hit_data Triangle::intersect(const Ray& ray) const {
  // May hold a_col1 and acol_2 as member?
  Hit_data hit_data;
  hit_data.t = std::numeric_limits<float>::infinity();
  hit_data.shape = NULL;
  /*if (!ray.shadow && ray.d.dot(normal) > 0.0f) {
    return hit_data;
  }*/
  const Vector3& v_0 = scene_->get_vertex_at(index_0);
  const Vector3& v_1 = scene_->get_vertex_at(index_1);
  const Vector3& v_2 = scene_->get_vertex_at(index_2);
  const Vector3 a_col1 = v_0 - v_1;
  const Vector3 a_col2 = v_0 - v_2;
  const Vector3& a_col3 = ray.d;
  const float det_a = determinant(a_col1, a_col2, a_col3);
  if (det_a > -kEpsilon && det_a < kEpsilon) {
    return hit_data;
  }
  const Vector3 b = (v_0 - ray.o) / det_a;
  const float beta = determinant(b, a_col2, a_col3);
  if (beta < 0.0f - kEpsilon || beta > 1.0f + kEpsilon) return hit_data;
  const float gamma = determinant(a_col1, b, a_col3);
  if (gamma < 0.0f - kEpsilon || beta + gamma > 1.0f + kEpsilon) {
    return hit_data;
  }
  const float t = determinant(a_col1, a_col2, b);
  if (t > -kEpsilon) {
    hit_data.t = t;
    hit_data.shape = this;
    hit_data.normal = this->normal;
  }
  return hit_data;
}
