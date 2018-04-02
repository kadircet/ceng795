#pragma once
#ifndef SPHERE_H_
#define SPHERE_H_
#include <cmath>
#include <limits>
#include "Ray.h"
#include "Shape.h"
#include "Vector3.h"
#include "Transformation.h"

class Sphere : public Shape {
 public:
  Vector3 center;
  float radius;
  int material_id;

  Sphere(const Vector3& center, float radius, int material_id, const Transformation& transformation)
      : center(center), radius(radius), material_id(material_id), transformation_(transformation) {
    is_identity_ = transformation_.get_transformation_matrix().is_identity();
    const Vector3 delta(radius);
    if (is_identity_) {
      bounding_box_ = Bounding_box(center - delta, center + delta);
    }
    else {
      bounding_box_ = Bounding_box::apply_transform(
        Bounding_box(center - delta, center + delta), 
        transformation_);
    }
   
  }

  int get_material_id() const override { return material_id; }
  const Bounding_box& get_bounding_box() const override {
    return bounding_box_;
  }
  bool intersect(const Ray& ray, Hit_data& hit_data) const override {
    if (is_identity_) {
      Vector3 center_to_origin = ray.o - center;
      const float a = ray.d.dot(ray.d);
      const float b = 2 * ray.d.dot(center_to_origin);
      const float c = center_to_origin.dot(center_to_origin) - radius * radius;
      const float determinant = b * b - 4 * a * c;
      if (determinant < -kEpsilon) {
        return false;
      }
      else if (determinant < kEpsilon) {
        hit_data.t = -b / (2 * a);
        hit_data.normal = (ray.point_at(hit_data.t) - center).normalize();
        hit_data.shape = this;
      }
      else {
        const float sqrt_det = sqrt(determinant);
        const float t1 = (-b + sqrt_det) / (2 * a);
        const float t2 = (-b - sqrt_det) / (2 * a);
        hit_data.t = fmin(t1, t2);
        if (t2 < .0) {
          hit_data.t = t1;
        }
        else {
          hit_data.t = t2;
        }
        hit_data.normal = (ray.point_at(hit_data.t) - center).normalize();
        hit_data.shape = this;
      }
      return true;
    }
    else {
      const Matrix4x4& inverse_transformation = transformation_.get_inverse_transformation_matrix();
      const Ray ray_local(inverse_transformation.multiply(ray.o),
        inverse_transformation.multiply(ray.d, true));
      Vector3 center_to_origin = ray_local.o - center;
      const float a = ray_local.d.dot(ray_local.d);
      const float b = 2 * ray_local.d.dot(center_to_origin);
      const float c = center_to_origin.dot(center_to_origin) - radius * radius;
      const float determinant = b * b - 4 * a * c;
      if (determinant < -kEpsilon) {
        return false;
      }
      else if (determinant < kEpsilon) {
        hit_data.t = -b / (2 * a);
        const Matrix4x4& normal_transformation = transformation_.get_normal_transformation_matrix();
        hit_data.normal = normal_transformation.multiply((ray_local.point_at(hit_data.t) - center).normalize(), true)
          .normalize();
        hit_data.shape = this;
      }
      else {
        const float sqrt_det = sqrt(determinant);
        const float t1 = (-b + sqrt_det) / (2 * a);
        const float t2 = (-b - sqrt_det) / (2 * a);
        hit_data.t = fmin(t1, t2);
        if (t2 < .0) {
          hit_data.t = t1;
        }
        else {
          hit_data.t = t2;
        }
        const Matrix4x4& normal_transformation = transformation_.get_normal_transformation_matrix();
        hit_data.normal = normal_transformation.multiply((ray_local.point_at(hit_data.t) - center).normalize(), true)
          .normalize();
        hit_data.shape = this;
      }
      return true;
    }
  }
  void print_debug(int indentation) const override {
    for (int index = 0; index < indentation; index++) {
      std::cout << "\t";
    }
    std::cout << "Sphere: " << center << "," << radius
              << " material: " << material_id << std::endl;
  }

 private:
  Transformation transformation_;
  bool is_identity_;
  Bounding_box bounding_box_;
};
#endif
