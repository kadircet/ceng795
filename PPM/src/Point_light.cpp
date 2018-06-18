#include "Point_light.h"
#include <sstream>
Point_light::Point_light(const Vector3& position, const Vector3& intensity)
    : position_(position), intensity_(intensity) {}

void Point_light::generate_photon(Ray& photon_ray, Vector3& flux,
                                  int photon_id) const {
  flux = intensity_ * (M_PI * 4.0f);
  float p = 2.0f * M_PI * hal(0, photon_id);
  float t = 2.0f * std::acos(clamp(0.0f, 1.0f, sqrt(1.0 - hal(1, photon_id))));
  float sint = sin(t);
  photon_ray.d = Vector3(std::cos(p) * sint, std::cos(t), std::sin(p) * sint);
  photon_ray.o = position_;
}

void Point_light::load_point_lights_from_xml(tinyxml2::XMLElement* element,
                                             std::vector<Light*>& lights) {
  element = element->FirstChildElement("PointLight");
  std::stringstream stream;
  while (element) {
    Vector3 position;
    Vector3 intensity;
    auto child = element->FirstChildElement("Position");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("Intensity");
    stream << child->GetText() << std::endl;

    stream >> position.x >> position.y >> position.z;
    stream >> intensity.x >> intensity.y >> intensity.z;

    lights.push_back(new Point_light(position, intensity));
    element = element->NextSiblingElement("PointLight");
  }
  stream.clear();
  std::cout << "PointLights are parsed" << std::endl;
}

/*Vector3 Point_light::direction_and_distance(const Vector3& from_point,
                                        const Vector3& normal,
                                        float& distance,
                                        float& probability) const {
const Vector3 direction = position_ - from_point;
distance = direction.length();
return direction;
}

Vector3 Point_light::incoming_radiance(const Vector3& from_point_to_light,
                                   float probability) const {
float x = from_point_to_light.x;
float y = from_point_to_light.y;
float z = from_point_to_light.z;
return intensity_ / (x * x + y * y + z * z);
}*/
