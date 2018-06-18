#include "Camera.h"
#include <sstream>
#include "Photographic_tmo.h"
void Camera::load_cameras_from_xml(tinyxml2::XMLElement* element,
                                   std::vector<Camera>& cameras) {
  constexpr float degrees_to_radians = M_PI / 180.0f;
  element = element->FirstChildElement("Camera");

  std::stringstream stream;
  while (element) {
    Vector3 position;
    Vector3 up;
    Vector3 gaze;
    float near_distance;
    float near_l, near_r, near_b, near_t;
    int image_width, image_height;
    int number_of_samples;
    std::string image_name;
    bool left_handed = false;

    auto child = element->FirstChildElement("Position");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("Up");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("NearDistance");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("ImageResolution");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("NumSamples");
    if (child) {
      stream << child->GetText() << std::endl;
    } else {
      stream << 1 << std::endl;
    }
    child = element->FirstChildElement("ImageName");
    stream << child->GetText() << std::endl;

    stream >> position.x >> position.y >> position.z;
    stream >> up.x >> up.y >> up.z;
    stream >> near_distance;
    stream >> image_width >> image_height;
    stream >> number_of_samples;
    number_of_samples = (int)sqrt(number_of_samples);
    if (number_of_samples <= 0) number_of_samples = 1;
    stream >> image_name;
    const char* camera_type = element->Attribute("type");
    if (camera_type && std::string(camera_type) == std::string("simple")) {
      child = element->FirstChildElement("Gaze");
      if (!child) {
        child = element->FirstChildElement("GazePoint");
      }
      stream << child->GetText() << std::endl;
      child = element->FirstChildElement("FovY");
      stream << child->GetText() << std::endl;

      Vector3 gaze_point;
      float FovY;
      stream >> gaze_point.x >> gaze_point.y >> gaze_point.z;
      stream >> FovY;
      float half_y_radian = degrees_to_radians * FovY / 2;
      near_t = tanf(half_y_radian) * near_distance;
      float aspect_ratio = 1.0f * image_width / image_height;
      near_b = -1.0f * near_t;
      near_r = near_t * aspect_ratio;
      near_l = -1.0f * near_r;
      gaze = (gaze_point - position).normalize();

    } else {
      child = element->FirstChildElement("Gaze");
      if (!child) {
        child = element->FirstChildElement("GazePoint");
      }
      stream << child->GetText() << std::endl;
      child = element->FirstChildElement("NearPlane");
      stream << child->GetText() << std::endl;
      stream >> gaze.x >> gaze.y >> gaze.z;
      stream >> near_l >> near_r >> near_b >> near_t;
    }
    Tonemapping_operator* tmo = nullptr;
    child = element->FirstChildElement("Tonemap");
    if (child) {
      auto tonemap_child = child->FirstChildElement("TMO");
      if (std::string(tonemap_child->GetText()) ==
          std::string("Photographic")) {
        float image_key, saturation_percentage, saturation;
        tonemap_child = child->FirstChildElement("TMOOptions");
        stream << tonemap_child->GetText() << std::endl;
        tonemap_child = child->FirstChildElement("Saturation");
        stream << tonemap_child->GetText() << std::endl;
        stream >> image_key >> saturation_percentage >> saturation;
        tmo =
            new Photographic_tmo(image_key, saturation_percentage, saturation);
      }
    }
    const char* handedness = element->Attribute("handedness");
    if (handedness && std::string(handedness) == std::string("left")) {
      left_handed = true;
    }
    cameras.push_back(
        std::move(Camera(up, gaze, position, number_of_samples, image_name,
                         near_l, near_r, near_b, near_t, near_distance,
                         image_width, image_height, tmo, left_handed)));
    element = element->NextSiblingElement("Camera");
  }
  stream.clear();
  std::cout << "Cameras are parsed" << std::endl;
}
