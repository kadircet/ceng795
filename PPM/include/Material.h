#pragma once
#ifndef MATERIAL_H_
#define MATERIAL_H_
#include <vector>
#include "Vector3.h"
#include "tinyxml2.h"
enum Material_type { mt_diffuse, mt_mirror, mt_refractive };
struct Material {
  static void load_materials_from_xml(tinyxml2::XMLElement* element,
                                      std::vector<Material>& materials);
  Vector3 diffuse;
  Vector3 specular;
  Vector3 mirror;
  Vector3 transparency;
  Material_type material_type;
  int brdf_id;
  float refraction_index;
  float phong_exponent;
};
#endif
