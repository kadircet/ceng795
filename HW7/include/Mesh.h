#pragma once
#ifndef MESH_H_
#define MESH_H_
#include <vector>
#include "Bounding_volume_hierarchy.h"
#include "Matrix4x4.h"
#include "Shape.h"
#include "Transformation.h"
#include "Triangle.h"

class Mesh : public Shape {
 public:
  int material_id;
  int texture_id;
  Shape* bvh;

  // Used for instances
  Vector3 velocity;
  Transformation base_transform;

  bool intersect(const Ray& ray, Hit_data& hit_data,
                 bool culling) const override {
    if (bvh->intersect(ray, hit_data, culling)) {
      hit_data.is_light_object = false;
      return true;
    }
    return false;
  }

  int get_material_id() const override { return material_id; }
  int get_texture_id() const override { return texture_id; }
  const Bounding_box& get_bounding_box() const override {
    return bvh->get_bounding_box();
  }

  Mesh(int material_id, int texture_id, std::vector<Shape*>& triangles,
       const Transformation& b_transform, const Vector3& velocity)
      : material_id(material_id),
        texture_id(texture_id),
        bvh(NULL),
        velocity(velocity),
        base_transform(b_transform) {
    bvh = BVH::create_bvh(triangles);
  }

  ~Mesh() {
    if (bvh) {
      delete bvh;
    }
  }
  void print_debug(int indentation) const override {
    for (int index = 0; index < indentation; index++) {
      std::cout << "\t";
    }
    std::cout << "Mesh->" << std::endl;
    bvh->print_debug(indentation);
  }
};

class Mesh_instance : public Shape {
 public:
  int material_id;
  int texture_id;
  Vector3 velocity;
  bool intersect(const Ray& ray, Hit_data& hit_data,
                 bool culling) const override {
    static const Vector3 zero_vector(0.0f);
    // Checking if ray hits the world space bounding box
    float bbox_t = bounding_box_.intersect(ray);
    if (bbox_t < 0.0f || bbox_t == kInf) {
      return false;
    }
    if (is_refractive_) {
      culling = false;
    }
    if (velocity == zero_vector) {
      const Matrix4x4& inverse_transformation =
          transformation_.get_inverse_transformation_matrix();
      Ray ray_local(inverse_transformation.multiply(ray.o),
                    inverse_transformation.multiply(ray.d, true), ray.ray_type,
                    ray.time);
      if (mesh_->intersect(ray_local, hit_data, culling)) {
        hit_data.normal = transformation_.get_normal_transformation_matrix()
                              .multiply(hit_data.normal, true)
                              .normalize();
        hit_data.is_light_object = false;
        hit_data.shape = this;
        return true;
      }
      return false;
    } else {
      // Motion blur
      Vector3 delta_position = ray.time * velocity;
      Translation translation(delta_position.x, delta_position.y,
                              delta_position.z);
      Matrix4x4 translated_transformation_matrix =
          translation.get_transformation_matrix() *
          transformation_.get_transformation_matrix();
      Arbitrary_transformation transformation(translated_transformation_matrix);
      const Matrix4x4& inverse_transformation =
          transformation.get_inverse_transformation_matrix();
      Ray ray_local(inverse_transformation.multiply(ray.o),
                    inverse_transformation.multiply(ray.d, true), ray.ray_type,
                    ray.time);
      if (mesh_->intersect(ray_local, hit_data, culling)) {
        hit_data.normal = transformation.get_normal_transformation_matrix()
                              .multiply(hit_data.normal, true)
                              .normalize();
        hit_data.shape = this;
        hit_data.is_light_object = false;
        return true;
      }
      return false;
    }
  }

  int get_material_id() const override { return material_id; }
  int get_texture_id() const override { return texture_id; }
  const Bounding_box& get_bounding_box() const override {
    return bounding_box_;
  }

  Mesh_instance(int material_id, int texture_id, const Mesh* mesh,
                const Transformation& transformation, const Vector3& velocity,
                bool is_refractive)
      : material_id(material_id),
        texture_id(texture_id),
        velocity(velocity),
        mesh_(mesh),
        transformation_(transformation),
        bounding_box_(Bounding_box::apply_transform(mesh->get_bounding_box(),
                                                    transformation)),
        is_refractive_(is_refractive) {
    if (velocity != Vector3(0.0f)) {
      Vector3 min = bounding_box_.min_corner;
      Vector3 max = bounding_box_.max_corner;
      bounding_box_.expand(Bounding_box(min + velocity, max + velocity));
      bounding_box_.expand(Bounding_box(min - velocity, max - velocity));
    }
  }

  void print_debug(int indentation) const override {}

 private:
  const Mesh* mesh_;
  const Transformation transformation_;
  Bounding_box bounding_box_;
  bool is_refractive_;
};
#endif
