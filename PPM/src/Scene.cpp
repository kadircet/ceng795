#include "Scene.h"
#include <cmath>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include "tinyxml2.h"
//#define GAUSSIAN_FILTER
Scene::Scene(const std::string& file_name) {
  tinyxml2::XMLDocument file;
  std::stringstream stream;
  auto res = file.LoadFile(file_name.c_str());
  if (res) {
    throw std::runtime_error("Error: The xml file cannot be loaded.");
  }
  auto root = file.FirstChild();
  if (!root) {
    throw std::runtime_error("Error: Root is not found.");
  }

  // Get ShadowRayEpsilon
  auto element = root->FirstChildElement("ShadowRayEpsilon");
  if (element) {
    stream << element->GetText() << std::endl;
  } else {
    stream << "0.001" << std::endl;
  }
  stream >> shadow_ray_epsilon;
  //

  // Get MaxRecursionDepth
  element = root->FirstChildElement("MaxRecursionDepth");
  if (element) {
    stream << element->GetText() << std::endl;
  } else {
    stream << "0" << std::endl;
  }
  stream >> max_recursion_depth;
  //

  // Get Cameras
  element = root->FirstChildElement("Cameras");
  if (element) {
    Camera::load_cameras_from_xml(element, cameras);
  }
  // Cameras End

  // Get Materials
  element = root->FirstChildElement("Materials");
  if (element) {
    Material::load_materials_from_xml(element, materials);
  }
  // Materials End

  // Get Transformations
  element = root->FirstChildElement("Transformations");
  if (element) {
    Translation::load_translation_transformations_from_xml(
        element, translation_transformations);
    Scaling::load_scaling_transformations_from_xml(element,
                                                   scaling_transformations);
    Rotation::load_rotation_transformations_from_xml(element,
                                                     rotation_transformations);
  }
  //
}

Scene::~Scene() {}
