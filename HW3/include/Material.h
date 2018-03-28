#pragma once
#ifndef MATERIAL_H_
#define MATERIAL_H_
#include "Vector3.h"
struct Material {
  Vector3 ambient;
  Vector3 diffuse;
  Vector3 specular;
  Vector3 mirror;
  Vector3 transparency;
  float refraction_index;
  float phong_exponent;
};
#endif
