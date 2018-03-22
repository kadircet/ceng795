#include "Scene.h"
#include <cmath>
#include <random>
#include <sstream>
#include <string>
#include "Pixel.h"
#include "tinyxml2.h"
inline int max(int a, int b) { return a > b ? a : b; }

inline int min(int a, int b) { return a < b ? a : b; }

inline float gaussian_filter(float x, float y, float sigma) {
  return exp(-(x * x + y * y) / (2 * sigma * sigma)) / (2 * M_PI * sigma);
}
void debug(const char* str) { std::cout << str << std::endl; }
void Scene::render_image(int camera_index, Pixel* result,
                         const int starting_row,
                         const int height_increase) const {
  const Camera& camera = cameras[camera_index];
  const Image_plane& image_plane = camera.get_image_plane();
  const int width = image_plane.width;
  const int height = image_plane.height;
  const int number_of_samples = camera.get_number_of_samples();
  if (number_of_samples == 1) {
    for (int j = starting_row; j < height; j += height_increase) {
      for (int i = 0; i < width; i++) {
        Vector3 color =
            trace_ray(camera.calculate_ray_at(i, j), max_recursion_depth);
        result[j * width + i].add_color(color, 1);
      }
    }
  } else {
    for (int j = starting_row; j < height; j += height_increase) {
      for (int i = 0; i < width; i++) {
        std::default_random_engine generator;
        generator.seed(
            std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_real_distribution<float> distribution(0.0, 1);
        for (int x = 0; x < number_of_samples; x++) {
          for (int y = 0; y < number_of_samples; y++) {
            float epsilon_x = distribution(generator);
            float epsilon_y = distribution(generator);
            float sample_x = (x + epsilon_x) / number_of_samples;
            float sample_y = (y + epsilon_y) / number_of_samples;
            // since calculate_ray_at adds 0.5 to the pixel number, subtracting
            // 0.5
            Vector3 color = trace_ray(
                camera.calculate_ray_at(i + sample_x - 0.5, j + sample_y - 0.5),
                max_recursion_depth);
            for (int affected_j = j - 1; affected_j < j + 2; affected_j++) {
              if (affected_j < 0 || affected_j >= height) {
                continue;
              }
              for (int affected_i = i - 1; affected_i < i + 2; affected_i++) {
                if (affected_i < 0 || affected_i >= width) {
                  continue;
                }
                float s_x = (i + sample_x) - (affected_i + 0.5f);
                float s_y = (j + sample_y) - (affected_j + 0.5f);
                result[affected_j * width + affected_i].add_color(
                    color, gaussian_filter(s_x, s_y, 1.0f / 3.0f));
              }
            }
          }
        }
      }
    }
  }
}
const Vector3 zero_vector(0.0f);
bool Scene::refract_ray(const Vector3& direction_unit, const Vector3& normal,
                        const float refraction_index,
                        Vector3& transmitted_d) const {
  float n_ratio = 1 / refraction_index;
  float cos_theta = (-direction_unit).dot(normal);
  float delta = 1 - (n_ratio) * (n_ratio) * (1 - (cos_theta * cos_theta));
  if (delta < 0.0f) {
    return false;
  }
  // TODO check if it is needed to be normalized.
  transmitted_d =
      ((direction_unit + normal * cos_theta) * n_ratio - normal * sqrt(delta))
          .normalize();
  return true;
}
Vector3 Scene::trace_ray(const Ray& ray, int current_recursion_depth) const {
  Vector3 color;
  // Find intersection
  Hit_data hit_data;
  hit_data.t = std::numeric_limits<float>::infinity();
  hit_data.shape = NULL;
  if (!bvh->intersect(ray, hit_data)) {
    // Check if it is primary ray or mirror ray
    if (this->max_recursion_depth == current_recursion_depth) {
      return background_color;
    }
    return color;
  }
  const Vector3 intersection_point = ray.point_at(hit_data.t);
  const Material& material = materials[hit_data.shape->get_material_id()];

  const Vector3 w_0 = (ray.o - intersection_point).normalize();
  const Vector3& normal = hit_data.normal;

  if (!ray.in_medium) {
    // ambient light
    // Should it be outside of if statement?
    color += material.ambient * ambient_light;

    // point lights
    for (const Point_light& point_light : point_lights) {
      // Shadow check
      const Vector3 light_distance_vec =
          point_light.position - intersection_point;
      const Vector3 w_i = light_distance_vec.normalize();
      float light_distance = light_distance_vec.length();
      Ray shadow_ray(intersection_point + (shadow_ray_epsilon * w_i), w_i);
      Hit_data shadow_hit_data;
      shadow_hit_data.t = std::numeric_limits<float>::infinity();
      shadow_hit_data.shape = NULL;
      bvh->intersect(shadow_ray, shadow_hit_data);
      if (shadow_hit_data.t < (light_distance - shadow_ray_epsilon) &&
          shadow_hit_data.t > 0.0f) {
        continue;
      }
      //
      float light_distance_squared = light_distance * light_distance;
      float diffuse_cos_theta = hit_data.normal.dot(w_i);
      color += material.diffuse * point_light.intensity * diffuse_cos_theta /
               light_distance_squared;
      float specular_cos_theta =
          fmax(0.0f, hit_data.normal.dot((w_0 + w_i).normalize()));
      color += material.specular * point_light.intensity *
               pow(specular_cos_theta, material.phong_exponent) /
               light_distance_squared;
    }
  }
  // Reflection
  if (material.mirror != zero_vector && current_recursion_depth > 0) {
    const Vector3 w_r = ((2 * normal.dot(w_0) * normal) - w_0).normalize();
    Ray mirror_ray(intersection_point + (w_r * shadow_ray_epsilon), w_r);
    color +=
        material.mirror * trace_ray(mirror_ray, current_recursion_depth - 1);
  }

  // Refraction
  if (material.transparency != zero_vector && current_recursion_depth > 0) {
    const Vector3 w_r = ((2 * normal.dot(w_0) * normal) - w_0).normalize();
    Vector3 transmission_direction = zero_vector;
    float cos_theta = 0.0f;
    Vector3 k(0.0f);
    Vector3 d_n = ray.d.normalize();
    float n = material.refraction_index;
    bool total_internal_reflection = false;
    bool entering_ray;
    if (d_n.dot(normal) < 0.0f && !ray.in_medium) {
      refract_ray(d_n, normal, n, transmission_direction);
      cos_theta = -d_n.dot(normal);
      k = Vector3(1.0f);
      entering_ray = true;
    } else {
      const Vector3& transparency = material.transparency;
      const float hit_data_t = hit_data.t;
      k.x = exp(-log(transparency.x) * hit_data_t);
      k.y = exp(-log(transparency.y) * hit_data_t);
      k.z = exp(-log(transparency.z) * hit_data_t);
      entering_ray = false;
      if (refract_ray(d_n, -normal, 1.0f / n, transmission_direction)) {
        cos_theta = transmission_direction.dot(normal);
      } else {
        total_internal_reflection = true;
      }
    }

    if (total_internal_reflection) {
      Ray reflection_ray(intersection_point + (w_r * shadow_ray_epsilon), w_r);
      reflection_ray.in_medium = true;
      color += k * trace_ray(reflection_ray, current_recursion_depth - 1);
    } else {
      float r_0 = (n - 1) * (n - 1) / ((n + 1) * (n + 1));
      float r = r_0 + (1 - r_0) * pow(1 - cos_theta, 5);
      Ray reflection_ray(intersection_point + (w_r * shadow_ray_epsilon), w_r);
      reflection_ray.in_medium = !entering_ray;
      Ray transmission_ray(
          intersection_point + (transmission_direction * shadow_ray_epsilon),
          transmission_direction);
      transmission_ray.in_medium = entering_ray;
      color += k * (r * trace_ray(reflection_ray, current_recursion_depth - 1) +
                    (1 - r) * trace_ray(transmission_ray,
                                        current_recursion_depth - 1));
    }
  }
  return color;
}

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
    child = element->FirstChildElement("NumSamples");
    if (child) {
      stream << child->GetText() << std::endl;
    } else {
      stream << 1 << std::endl;
    }
    child = element->FirstChildElement("ImageName");
    stream << child->GetText() << std::endl;
    Vector3 position, up, gaze;
    float near_distance;
    float near_l, near_r, near_b, near_t;
    int image_width, image_height;
    int number_of_samples;
    std::string image_name;
    stream >> position.x >> position.y >> position.z;
    stream >> gaze.x >> gaze.y >> gaze.z;
    stream >> up.x >> up.y >> up.z;
    stream >> near_l >> near_r >> near_b >> near_t;
    stream >> near_distance;
    stream >> image_width >> image_height;
    stream >> number_of_samples;
    number_of_samples = sqrt(number_of_samples);
    if (number_of_samples <= 0) number_of_samples = 1;
    stream >> image_name;
    Camera camera(up, gaze, position, number_of_samples, image_name, near_l,
                  near_r, near_b, near_t, near_distance, image_width,
                  image_height);
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
    if (child) {
      stream << child->GetText() << std::endl;
    } else {
      stream << "0 0 0" << std::endl;
    }
    child = element->FirstChildElement("DiffuseReflectance");
    if (child) {
      stream << child->GetText() << std::endl;
    } else {
      stream << "0 0 0" << std::endl;
    }
    child = element->FirstChildElement("SpecularReflectance");
    if (child) {
      stream << child->GetText() << std::endl;
    } else {
      stream << "0 0 0" << std::endl;
    }
    child = element->FirstChildElement("MirrorReflectance");
    if (child) {
      stream << child->GetText() << std::endl;
    } else {
      stream << "0 0 0" << std::endl;
    }
    child = element->FirstChildElement("PhongExponent");
    if (child) {
      stream << child->GetText() << std::endl;
    } else {
      stream << "1" << std::endl;
    }
    child = element->FirstChildElement("Transparency");
    if (child) {
      stream << child->GetText() << std::endl;
    } else {
      stream << "0 0 0" << std::endl;
    }
    child = element->FirstChildElement("RefractionIndex");
    if (child) {
      stream << child->GetText() << std::endl;
    } else {
      stream << "1.0" << std::endl;
    }

    stream >> material.ambient.x >> material.ambient.y >> material.ambient.z;
    stream >> material.diffuse.x >> material.diffuse.y >> material.diffuse.z;
    stream >> material.specular.x >> material.specular.y >> material.specular.z;
    stream >> material.mirror.x >> material.mirror.y >> material.mirror.z;
    stream >> material.phong_exponent;
    stream >> material.transparency.x >> material.transparency.y >>
        material.transparency.z;
    stream >> material.refraction_index;

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

  std::vector<Shape*> objects;
  // Get Meshes
  element = root->FirstChildElement("Objects");
  element = element->FirstChildElement("Mesh");

  while (element) {
    child = element->FirstChildElement("Material");
    int material_id;
    stream << child->GetText() << std::endl;
    stream >> material_id;
    material_id--;

    child = element->FirstChildElement("Faces");
    stream << child->GetText() << std::endl;
    int v0_id, v1_id, v2_id;
    std::vector<Shape*> triangles;
    while (!(stream >> v0_id).eof()) {
      stream >> v1_id >> v2_id;
      triangles.push_back(
          new Triangle(this, v0_id - 1, v1_id - 1, v2_id - 1, material_id));
    }
    stream.clear();

    objects.push_back(new Mesh(material_id, -1, triangles));
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
    objects.push_back(
        new Triangle(this, v0_id - 1, v1_id - 1, v2_id - 1, material_id));
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

    objects.push_back(new Sphere(center_of_sphere, radius, material_id));
    element = element->NextSiblingElement("Sphere");
  }
  bvh = BVH::create_bvh(objects);
  // bvh->print_debug(0);
}

Scene::~Scene() { delete bvh; }
