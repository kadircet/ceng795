#include "Triangle.h"
#include <algorithm>
#include "Hit_data.h"
#include "Scene.h"
Triangle::Triangle(const Scene* scene, int index_0, int index_1, int index_2,
                   int offset, int material_id, int texture_id,
                   const Transformation& transformation)
    : index_0(index_0),
      index_1(index_1),
      index_2(index_2),
      offset(offset),
      material_id(material_id),
      texture_id(texture_id),
      transformation_(transformation),
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
  is_identity_ = transformation_.get_transformation_matrix().is_identity();
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
  if (is_identity_) {
    bounding_box_ = Bounding_box(min_c, max_c);
  } else {
    bounding_box_ = Bounding_box::apply_transform(Bounding_box(min_c, max_c),
                                                  transformation_);
  }
}
float Triangle::get_surface_area() const {
  const Vector3& v_0 =
      scene_->get_vertex_at(index_0 + offset).get_vertex_position();
  const Vector3& v_1 =
      scene_->get_vertex_at(index_1 + offset).get_vertex_position();
  const Vector3& v_2 =
      scene_->get_vertex_at(index_2 + offset).get_vertex_position();
  return (v_1 - v_0).cross(v_2 - v_0).length() / 2;
}
bool Triangle::intersect(const Ray& ray, Hit_data& hit_data) const {
  // May hold a_col1 and acol_2 as member?
  /*if (!ray.shadow && ray.d.dot(normal) > 0.0f) {
    return hit_data;
  }*/
  const Vector3& v_0 =
      scene_->get_vertex_at(index_0 + offset).get_vertex_position();
  const Vector3& v_1 =
      scene_->get_vertex_at(index_1 + offset).get_vertex_position();
  const Vector3& v_2 =
      scene_->get_vertex_at(index_2 + offset).get_vertex_position();
  const Vector3 a_col1 = v_0 - v_1;
  const Vector3 a_col2 = v_0 - v_2;
  if (is_identity_) {
    const Vector3& a_col3 = ray.d;
    const float det_a = determinant(a_col1, a_col2, a_col3);
    if (det_a == 0.0f) {
      return false;
    }
    const Vector3 b = (v_0 - ray.o) / det_a;
    const float beta = determinant(b, a_col2, a_col3);
    if (beta < 0.0f || beta > 1.0f) return false;
    const float gamma = determinant(a_col1, b, a_col3);
    if (gamma < 0.0f || beta + gamma > 1.0f) {
      return false;
    }
    const float t = determinant(a_col1, a_col2, b);
    if (t > 0.0f) {
      hit_data.t = t;
      hit_data.shape = this;
      hit_data.normal = this->normal;
      return true;
    }
    return false;
  } else {
    const Matrix4x4& inverse_transformation =
        transformation_.get_inverse_transformation_matrix();
    const Ray ray_local(inverse_transformation.multiply(ray.o),
                        inverse_transformation.multiply(ray.d, true),
                        ray.ray_type);
    const Vector3& a_col3 = ray_local.d;
    const float det_a = determinant(a_col1, a_col2, a_col3);
    if (det_a == 0.0f) {
      return false;
    }
    const Vector3 b = (v_0 - ray_local.o) / det_a;
    const float beta = determinant(b, a_col2, a_col3);
    if (beta < 0.0f || beta > 1.0f) return false;
    const float gamma = determinant(a_col1, b, a_col3);
    if (gamma < 0.0f || beta + gamma > 1.0f) {
      return false;
    }
    const float t = determinant(a_col1, a_col2, b);
    if (t > 0.0f) {
      hit_data.t = t;
      hit_data.shape = this;
      const Matrix4x4& normal_transformation =
          transformation_.get_normal_transformation_matrix();
      // TODO: Check if it is precomputable?
      hit_data.normal =
          normal_transformation.multiply(this->normal, true).normalize();
      return true;
    }
    // TODO ADD Texture support for primitive triangles
    return false;
  }
}
