#pragma once

#ifndef TRANSFORMATION_H_
#define TRANSFORMATION_H_
#include <vector>
#include "Matrix4x4.h"
#include "tinyxml2.h"

class Transformation {
 public:
  const Matrix4x4& get_transformation_matrix() const { return transformation_; }
  const Matrix4x4& get_normal_transformation_matrix() const {
    return normal_transformation_;
  }
  const Matrix4x4& get_inverse_transformation_matrix() const {
    return inverse_transformation_;
  }

 protected:
  Matrix4x4 transformation_;
  Matrix4x4 normal_transformation_;
  Matrix4x4 inverse_transformation_;
};

class Scaling : public Transformation {
 public:
  Scaling(float x, float y, float z);
  static void load_scaling_transformations_from_xml(
      tinyxml2::XMLElement* element,
      std::vector<Scaling>& scaling_transformations);
};

class Translation : public Transformation {
 public:
  Translation(float x, float y, float z);
  static void load_translation_transformations_from_xml(
      tinyxml2::XMLElement* element,
      std::vector<Translation>& translation_transformations);
};

class Rotation : public Transformation {
 public:
  Rotation(float angle, float x, float y, float z);
  static void load_rotation_transformations_from_xml(
      tinyxml2::XMLElement* element,
      std::vector<Rotation>& rotation_transformations);
};

class Arbitrary_transformation : public Transformation {
 public:
  Arbitrary_transformation(const Matrix4x4& transformation_matrix);
};

#endif
