#include "Scene.h"
#include <cmath>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include "Bounding_volume_hierarchy.h"
#include "Mesh.h"
#include "Sphere.h"
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
  std::cout << "ShadowRayEpsilon is parsed" << std::endl;
  //

  // Get MaxRecursionDepth
  element = root->FirstChildElement("MaxRecursionDepth");
  if (element) {
    stream << element->GetText() << std::endl;
  } else {
    stream << "0" << std::endl;
  }
  stream >> max_recursion_depth;
  std::cout << "MaxRecursionDepth is parsed" << std::endl;
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
  // Transformations End

  // Get VertexData
  element = root->FirstChildElement("VertexData");
  if (element) {
    const char* binary_file = element->Attribute("binaryFile");
    if (binary_file) {
      // parse_binary_vertexdata(std::string(binary_file));
    } else {
      stream << element->GetText() << std::endl;
      Vector3 vertex;
      while (!(stream >> vertex.x).eof()) {
        stream >> vertex.y >> vertex.z;
        vertex_data.push_back(Vertex(vertex));
      }
    }
  }
  stream.clear();
  std::cout << "VertexData is parsed" << std::endl;
  // VertexData End

  // Get Objects
  std::vector<Shape*> objects;
  element = root->FirstChildElement("Objects");
  if (element) {
    Sphere::load_spheres_from_xml(
        this, element, objects, vertex_data, scaling_transformations,
        translation_transformations, rotation_transformations);
    Mesh::load_meshes_from_xml(
        this, element, meshes, vertex_data, scaling_transformations,
        translation_transformations, rotation_transformations);
    Mesh_instance::create_mesh_instances_for_meshes(meshes, objects, materials);
    Mesh_instance::load_mesh_instances_from_xml(
        element, meshes, objects, materials, scaling_transformations,
        translation_transformations, rotation_transformations);
  }
  for (Vertex& vertex : vertex_data) {
    vertex.finalize_normal();
  }
  bvh = BVH::create_bvh(objects);

  //
}

Scene::~Scene() {}
