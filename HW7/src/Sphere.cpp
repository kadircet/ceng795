#include "Sphere.h"
#include <algorithm>
#include "Ray.h"
#include "Scene.h"
#include "Texture.h"
Sphere::Sphere(const Scene* scene, const Vector3& center, float radius,
               int material_id, int texture_id,
               const Transformation& transformation, const Vector3& velocity)
    : center(center),
      radius(radius),
      material_id(material_id),
      texture_id(texture_id),
      velocity(velocity),
      transformation_(transformation),
      scene_(scene) {
  is_identity_ = transformation_.get_transformation_matrix().is_identity();
  const Vector3 delta(radius);
  if (is_identity_) {
    bounding_box_ = Bounding_box(center - delta, center + delta);
  } else {
    bounding_box_ = Bounding_box::apply_transform(
        Bounding_box(center - delta, center + delta), transformation_);
  }
  if (velocity != Vector3(0.0f)) {
    Vector3 min = bounding_box_.min_corner;
    Vector3 max = bounding_box_.max_corner;
    bounding_box_.expand(Bounding_box(min + velocity, max + velocity));
    bounding_box_.expand(Bounding_box(min - velocity, max - velocity));
  }
}

bool Sphere::intersect(const Ray& ray, Hit_data& hit_data) const {
  static const Vector3 zero_vector(0.0f);
  Vector3 delta_position = ray.time * velocity;
  Translation translation(delta_position.x, delta_position.y, delta_position.z);
  Matrix4x4 translated_transformation_matrix =
      translation.get_transformation_matrix() *
      transformation_.get_transformation_matrix();
  Arbitrary_transformation transformation(translated_transformation_matrix);
  const Matrix4x4& inverse_transformation =
      transformation.get_inverse_transformation_matrix();
  const Ray ray_local(inverse_transformation.multiply(ray.o),
                      inverse_transformation.multiply(ray.d, true),
                      ray.ray_type);
  Vector3 center_to_origin = ray_local.o - center;
  const float a = ray_local.d.dot(ray_local.d);
  const float b = 2 * ray_local.d.dot(center_to_origin);
  const float c = center_to_origin.dot(center_to_origin) - radius * radius;
  const float determinant = b * b - 4 * a * c;
  if (determinant < -intersection_test_epsilon) {
    return false;
  } else if (determinant < intersection_test_epsilon) {
    hit_data.t = -b / (2 * a);
  } else {
    const float sqrt_det = sqrt(determinant);
    const float t1 = (-b + sqrt_det) / (2 * a);
    const float t2 = (-b - sqrt_det) / (2 * a);
    if (t2 < 0.0f) {
      hit_data.t = t1;
    } else {
      hit_data.t = t2;
    }
  }
  const Matrix4x4& normal_transformation =
      transformation.get_normal_transformation_matrix();
  Vector3 local_intersection_point = ray_local.point_at(hit_data.t);
  Vector3 local_coordinates = local_intersection_point - center;
  Vector3 normal = local_coordinates.normalize();
  float u = -1;
  float v = -1;
  float perlin_value = -1;
  if (texture_id != -1) {
    get_uv(local_coordinates, u, v);
    const Texture& texture = scene_->textures[texture_id];
    if (texture.is_perlin_noise()) {
      perlin_value =
          texture.get_perlin_noise()->get_value_at(local_intersection_point);
      if (texture.is_bump()) {
        normal = texture.bump_normal(normal, local_intersection_point);
      }
    } else {
      if (texture.is_bump()) {
        // calculate gradient vectors
        float theta = M_PI * v;
        float phi = M_PI * (1 - 2 * u);
        Vector3 dp_du(local_coordinates.z * 2 * M_PI, 0,
                      local_coordinates.x * -2 * M_PI);
        Vector3 dp_dv(local_coordinates.y * M_PI * cos(phi),
                      -radius * M_PI * sin(theta),
                      local_coordinates.y * M_PI * sin(phi));
        normal = texture.bump_normal(dp_du, dp_dv, normal, u, v);
      }
    }
  }
  hit_data.u = u;
  hit_data.v = v;
  hit_data.perlin_value = perlin_value;
  hit_data.normal = normal_transformation.multiply(normal, true).normalize();

  hit_data.shape = this;
  hit_data.is_light_object = false;
  return true;
}

void Sphere::get_uv(const Vector3& local_coordinates, float& u,
                    float& v) const {
  float theta =
      std::acos(std::max(-1.0f, std::min(1.0f, local_coordinates.y / radius)));
  float phi = atan2(local_coordinates.z, local_coordinates.x);
  u = (M_PI - phi) / (2 * M_PI);
  v = theta / M_PI;
}
