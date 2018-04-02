#pragma once
#ifndef MESH_H_
#define MESH_H_
#include <vector>
#include "Bounding_volume_hierarchy.h"
#include "Shape.h"
#include "Triangle.h"
#include "Transformation.h"
#include "Matrix4x4.h"

class Mesh : public Shape {
 public:
  int material_id;
  int texture_id;
  Shape* bvh;

  //Used for instances
  Transformation base_transform;


  bool intersect(const Ray& ray, Hit_data& hit_data) const override {
    return bvh->intersect(ray, hit_data);
  }

  int get_material_id() const override { return material_id; }
  const Bounding_box& get_bounding_box() const override {
    return bvh->get_bounding_box();
  }

  Mesh(int material_id, int texture_id, std::vector<Shape*>& triangles, const Transformation& b_transform)
      : material_id(material_id), texture_id(texture_id), bvh(NULL), base_transform(b_transform) {
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
	bool intersect(const Ray& ray, Hit_data& hit_data) const override {
		//Checking if ray hits the world space bounding box
		float bbox_t = bounding_box_.intersect(ray);
		if (bbox_t < 0.0f || bbox_t == kInf) {
			return false;
		}

		const Matrix4x4& inverse_tranformation = transformation_.get_inverse_transformation_matrix();
		Ray ray_local(inverse_tranformation.multiply(ray.o),
			inverse_tranformation.multiply(ray.d, true));
		if (mesh_->intersect(ray_local, hit_data))
		{
			hit_data.normal = transformation_.get_normal_transformation_matrix().multiply(hit_data.normal, true).normalize();
      hit_data.shape = this;
			return true;
		}
		return false;
	}

	int get_material_id() const override { return material_id; }
	const Bounding_box& get_bounding_box() const override {
		return bounding_box_;
	}

	Mesh_instance(int material_id, int texture_id, const Mesh* mesh, const Transformation& transformation)
		: material_id(material_id), texture_id(texture_id), mesh_(mesh), transformation_(transformation),
		bounding_box_(Bounding_box::apply_transform(mesh->get_bounding_box(), transformation))
	{
	}
	
	void print_debug(int indentation) const override {

	}
private:
	const Mesh* mesh_;
	const Transformation transformation_;
	const Bounding_box bounding_box_;
};
#endif
