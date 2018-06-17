#include "Sphere.h"
#include <algorithm>
#include <sstream>
#include "Ray.h"
#include "Scene.h"
//#include "Texture.h"
Sphere::Sphere(const Scene* scene, const Vector3& center, float radius,
               int material_id, int texture_id,
               const Transformation& transformation)
    : center(center),
      radius(radius),
      material_id(material_id),
      texture_id(texture_id),
      transformation_(transformation),
      scene_(scene) {
  bool is_identity = transformation_.get_transformation_matrix().is_identity();
  const Vector3 delta(radius);
  if (is_identity) {
    bounding_box_ = Bounding_box(center - delta, center + delta);
  } else {
    bounding_box_ = Bounding_box::apply_transform(
        Bounding_box(center - delta, center + delta), transformation_);
  }
}

bool Sphere::intersect(const Ray& ray, Intersection& intersection,
                       bool culling) const {
  const Matrix4x4& inverse_transformation =
      transformation_.get_inverse_transformation_matrix();
  const Ray ray_local(inverse_transformation.multiply(ray.o),
                      inverse_transformation.multiply(ray.d, true),
                      ray.ray_type);
  Vector3 center_to_origin = ray_local.o - center;
  const float a = ray_local.d.dot(ray_local.d);
  const float b = 2 * ray_local.d.dot(center_to_origin);
  const float c = center_to_origin.dot(center_to_origin) - radius * radius;
  const float determinant = b * b - 4 * a * c;
  if (determinant < -intersection_test_epsilon) {
    return false;
  } else if (determinant < intersection_test_epsilon) {
    intersection.t = -b / (2 * a);
  } else {
    const float sqrt_det = sqrt(determinant);
    const float t1 = (-b + sqrt_det) / (2 * a);
    const float t2 = (-b - sqrt_det) / (2 * a);
    if (t2 < 0.0f) {
      intersection.t = t1;
    } else {
      intersection.t = t2;
    }
  }
  const Matrix4x4& normal_transformation =
      transformation_.get_normal_transformation_matrix();
  Vector3 local_intersection_point = ray_local.point_at(intersection.t);
  Vector3 local_coordinates = local_intersection_point - center;
  Vector3 normal = local_coordinates.normalize();
  intersection.normal =
      normal_transformation.multiply(normal, true).normalize();
  intersection.shape = this;
  return true;
}

void Sphere::get_uv(const Vector3& local_coordinates, float& u,
                    float& v) const {
  float theta =
      std::acos(std::max(-1.0f, std::min(1.0f, local_coordinates.y / radius)));
  float phi = atan2(local_coordinates.z, local_coordinates.x);
  u = (M_PI - phi) / (2 * M_PI);
  v = theta / M_PI;
}

void Sphere::load_spheres_from_xml(
    const Scene* scene, tinyxml2::XMLElement* element,
    std::vector<Shape*>& objects, std::vector<Vertex>& vertex_data,
    std::vector<Scaling>& scaling_transformations,
    std::vector<Translation>& translation_transformations,
    std::vector<Rotation>& rotation_transformations) {
  element = element->FirstChildElement("Sphere");
  std::stringstream stream;
  while (element) {
    int material_id;
    auto child = element->FirstChildElement("Material");
    stream << child->GetText() << std::endl;
    stream >> material_id;
    material_id--;

    child = element->FirstChildElement("Center");
    stream << child->GetText() << std::endl;
    int center;
    stream >> center;
    const Vector3& center_of_sphere =
        vertex_data[center - 1].get_vertex_position();

    float radius;
    child = element->FirstChildElement("Radius");
    stream << child->GetText() << std::endl;
    stream >> radius;
    Matrix4x4 arbitrary_transformation(true);
    child = element->FirstChildElement("Transformations");
    if (child) {
      char type;
      int index;
      stream.clear();
      stream << child->GetText() << std::endl;
      while (!(stream >> type).eof()) {
        stream >> index;
        index--;
        switch (type) {
          case 's':
            arbitrary_transformation =
                scaling_transformations[index].get_transformation_matrix() *
                arbitrary_transformation;
            break;
          case 't':
            arbitrary_transformation =
                translation_transformations[index].get_transformation_matrix() *
                arbitrary_transformation;
            break;
          case 'r':
            arbitrary_transformation =
                rotation_transformations[index].get_transformation_matrix() *
                arbitrary_transformation;
            break;
        }
      }
      stream.clear();
    }
    stream.clear();
    int texture_id = -1;
    child = element->FirstChildElement("Texture");
    if (child) {
      stream << child->GetText() << std::endl;
      stream >> texture_id;
      texture_id--;
    }
    stream.clear();
    objects.push_back(
        new Sphere(scene, center_of_sphere, radius, material_id, texture_id,
                   Arbitrary_transformation(arbitrary_transformation)));
    element = element->NextSiblingElement("Sphere");
  }
  stream.clear();
  std::cout << "Spheres are parsed" << std::endl;
}
