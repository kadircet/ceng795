#pragma once
#ifndef MESH_H_
#define MESH_H_
#include <vector>
#include "Bounding_volume_hierarchy.h"
#include "Material.h"
#include "Matrix4x4.h"
#include "Shape.h"
#include "Transformation.h"
#include "Vertex.h"
#include "tinyxml2.h"
class Scene;
class Mesh : public Shape {
 public:
  int material_id;
  int texture_id;
  Shape* bvh;

  // Used for creating instances
  Transformation base_transform;

  bool intersect(const Ray& ray, Intersection& intersection,
                 bool culling) const override {
    if (bvh->intersect(ray, intersection, culling)) {
      // intersection.is_light_object = false;
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
       const Transformation& b_transform)
      : material_id(material_id),
        texture_id(texture_id),
        bvh(NULL),
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
  static void load_meshes_from_xml(
      const Scene* scene, tinyxml2::XMLElement* element,
      std::vector<Mesh*>& meshes, std::vector<Vertex>& vertex_data,
      std::vector<Scaling>& scaling_transformations,
      std::vector<Translation>& translation_transformations,
      std::vector<Rotation>& rotation_transformations);
};

class Mesh_instance : public Shape {
 public:
  int material_id;
  int texture_id;
  bool intersect(const Ray& ray, Intersection& intersection,
                 bool culling) const override {
    // Checking if ray hits the world space bounding box
    float bbox_t = bounding_box_.intersect(ray);
    if (bbox_t < 0.0f || bbox_t == kInf) {
      return false;
    }
    if (is_refractive_) {
      culling = false;
    }

    const Matrix4x4& inverse_transformation =
        transformation_.get_inverse_transformation_matrix();
    Ray ray_local(inverse_transformation.multiply(ray.o),
                  inverse_transformation.multiply(ray.d, true), ray.ray_type);
    if (mesh_->intersect(ray_local, intersection, culling)) {
      intersection.normal = transformation_.get_normal_transformation_matrix()
                                .multiply(intersection.normal, true)
                                .normalize();
      // intersection.is_light_object = false;
      intersection.shape = this;
      return true;
    }
    return false;
  }

  int get_material_id() const override { return material_id; }
  int get_texture_id() const override { return texture_id; }
  const Bounding_box& get_bounding_box() const override {
    return bounding_box_;
  }

  Mesh_instance(int material_id, int texture_id, const Mesh* mesh,
                const Transformation& transformation, bool is_refractive)
      : material_id(material_id),
        texture_id(texture_id),
        mesh_(mesh),
        transformation_(transformation),
        bounding_box_(Bounding_box::apply_transform(mesh->get_bounding_box(),
                                                    transformation)),
        is_refractive_(is_refractive) {}

  void print_debug(int indentation) const override {}
  static void create_mesh_instances_for_meshes(
      std::vector<Mesh*>& meshes, std::vector<Shape*>& objects,
      std::vector<Material>& materials);
  static void load_mesh_instances_from_xml(
      tinyxml2::XMLElement* element, std::vector<Mesh*>& meshes,
      std::vector<Shape*>& objects, std::vector<Material>& materials,
      std::vector<Scaling>& scaling_transformations,
      std::vector<Translation>& translation_transformations,
      std::vector<Rotation>& rotation_transformations);

 private:
  const Mesh* mesh_;
  const Transformation transformation_;
  Bounding_box bounding_box_;
  bool is_refractive_;
};
#endif
