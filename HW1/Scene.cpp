#include "Scene.h"
#include <sstream>
#include <string>
#include "tinyxml2.h"
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

  // Get BackgroundColor
  auto element = root->FirstChildElement("BackgroundColor");
  if (element) {
    stream << element->GetText() << std::endl;
  } else {
    stream << "0 0 0" << std::endl;
  }
  stream >> background_color.x >> background_color.y >> background_color.z;

  // Get ShadowRayEpsilon
  element = root->FirstChildElement("ShadowRayEpsilon");
  if (element) {
    stream << element->GetText() << std::endl;
  } else {
    stream << "0.001" << std::endl;
  }
  stream >> shadow_ray_epsilon;

  // Get MaxRecursionDepth
  element = root->FirstChildElement("MaxRecursionDepth");
  if (element) {
    stream << element->GetText() << std::endl;
  } else {
    stream << "0" << std::endl;
  }
  stream >> max_recursion_depth;

  // Get Cameras
  element = root->FirstChildElement("Cameras");
  element = element->FirstChildElement("Camera");
  while (element) {
    auto child = element->FirstChildElement("Position");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("Gaze");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("Up");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("NearPlane");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("NearDistance");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("ImageResolution");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("ImageName");
    stream << child->GetText() << std::endl;
    Vector3 position, up, gaze;
    float near_distance;
    float near_l, near_r, near_b, near_t;
    int image_width, image_height;
    std::string image_name;
    stream >> position.x >> position.y >> position.z;
    stream >> gaze.x >> gaze.y >> gaze.z;
    stream >> up.x >> up.y >> up.z;
    stream >> near_l >> near_r >> near_b >> near_t;
    stream >> near_distance;
    stream >> image_width >> image_height;
    stream >> image_name;
    Camera camera(up, gaze, position, image_name, near_l, near_r, near_b,
                  near_t, near_distance, image_width, image_height);
    cameras.push_back(camera);
    element = element->NextSiblingElement("Camera");
  }

  // Get Lights
  element = root->FirstChildElement("Lights");
  auto child = element->FirstChildElement("AmbientLight");
  stream << child->GetText() << std::endl;
  stream >> ambient_light.x >> ambient_light.y >> ambient_light.z;
  element = element->FirstChildElement("PointLight");
  Point_light point_light;
  while (element) {
    child = element->FirstChildElement("Position");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("Intensity");
    stream << child->GetText() << std::endl;

    stream >> point_light.position.x >> point_light.position.y >>
        point_light.position.z;
    stream >> point_light.intensity.x >> point_light.intensity.y >>
        point_light.intensity.z;

    point_lights.push_back(point_light);
    element = element->NextSiblingElement("PointLight");
  }

  // Get Materials
  element = root->FirstChildElement("Materials");
  element = element->FirstChildElement("Material");
  Material material;
  while (element) {
    child = element->FirstChildElement("AmbientReflectance");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("DiffuseReflectance");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("SpecularReflectance");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("MirrorReflectance");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("PhongExponent");
    stream << child->GetText() << std::endl;

    stream >> material.ambient.x >> material.ambient.y >> material.ambient.z;
    stream >> material.diffuse.x >> material.diffuse.y >> material.diffuse.z;
    stream >> material.specular.x >> material.specular.y >> material.specular.z;
    stream >> material.mirror.x >> material.mirror.y >> material.mirror.z;
    stream >> material.phong_exponent;

    materials.push_back(material);
    element = element->NextSiblingElement("Material");
  }

  // Get VertexData
  element = root->FirstChildElement("VertexData");
  stream << element->GetText() << std::endl;
  Vector3 vertex;
  while (!(stream >> vertex.x).eof()) {
    stream >> vertex.y >> vertex.z;
    vertex_data.push_back(vertex);
  }
  stream.clear();

  // Get Meshes
  element = root->FirstChildElement("Objects");
  element = element->FirstChildElement("Mesh");
  Mesh mesh;
  while (element) {
    child = element->FirstChildElement("Material");
    stream << child->GetText() << std::endl;
    stream >> mesh.material_id;
    mesh.material_id--;

    child = element->FirstChildElement("Faces");
    stream << child->GetText() << std::endl;
    int v0_id, v1_id, v2_id;
    while (!(stream >> v0_id).eof()) {
      stream >> v1_id >> v2_id;
      const Vector3& v_0 = vertex_data[v0_id - 1];
      const Vector3& v_1 = vertex_data[v1_id - 1];
      const Vector3& v_2 = vertex_data[v2_id - 1];
      mesh.triangles.push_back(Triangle(v_0, v_1, v_2, mesh.material_id));
    }
    stream.clear();

    meshes.push_back(std::move(mesh));
    mesh.triangles.clear();
    element = element->NextSiblingElement("Mesh");
  }
  stream.clear();

  // Get Triangles
  element = root->FirstChildElement("Objects");
  element = element->FirstChildElement("Triangle");
  while (element) {
    int material_id;
    child = element->FirstChildElement("Material");
    stream << child->GetText() << std::endl;
    stream >> material_id;
    material_id--;

    child = element->FirstChildElement("Indices");
    stream << child->GetText() << std::endl;
    int v0_id, v1_id, v2_id;
    stream >> v0_id >> v1_id >> v2_id;
    const Vector3& v_0 = vertex_data[v0_id - 1];
    const Vector3& v_1 = vertex_data[v1_id - 1];
    const Vector3& v_2 = vertex_data[v2_id - 1];
    triangles.push_back(Triangle(v_0, v_1, v_2, material_id));
    element = element->NextSiblingElement("Triangle");
  }

  // Get Spheres
  element = root->FirstChildElement("Objects");
  element = element->FirstChildElement("Sphere");
  while (element) {
    int material_id;
    child = element->FirstChildElement("Material");
    stream << child->GetText() << std::endl;
    stream >> material_id;
    material_id--;

    child = element->FirstChildElement("Center");
    stream << child->GetText() << std::endl;
    int center;
    stream >> center;
    const Vector3& center_of_sphere = vertex_data[center - 1];

    float radius;
    child = element->FirstChildElement("Radius");
    stream << child->GetText() << std::endl;
    stream >> radius;

    spheres.push_back(Sphere(center_of_sphere, radius, material_id));
    element = element->NextSiblingElement("Sphere");
  }
}
