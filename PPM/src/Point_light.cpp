#include "Point_light.h"
#include <random>
#include <sstream>
Point_light::Point_light(const Vector3& position, const Vector3& intensity)
    : position_(position), intensity_(intensity) {}

void Point_light::generate_photon(Ray& photon_ray, Vector3& flux) const {
  flux = intensity_ * (M_PI * 4.0f);
  // Sample sphere
  thread_local static std::random_device rd;
  thread_local static std::mt19937 generator(rd());
  std::uniform_real_distribution<float> uniform_dist(0.0f, 1.0f);
  float epsilon_1 = uniform_dist(generator);
  float epsilon_2 = uniform_dist(generator);

  float phi;
  float theta;
  Vector3 w(0.0f, 1.0f, 0.0f);
  const Vector3 u = ((w.x != 0.0f || w.y != 0.0f) ? Vector3(-w.y, w.x, 0.0f)
                                                  : Vector3(0.0f, 1.0f, 0.0f))
                        .normalize();
  const Vector3 v = w.cross(u);
  phi = 2 * M_PI * epsilon_1;
  theta = 2 * M_PI * epsilon_2;

  photon_ray.d = (w * std::cos(theta) + v * std::sin(theta) * std::cos(phi) +
                  u * std::sin(theta) * std::sin(phi))
                     .normalize();
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
