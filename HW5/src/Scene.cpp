#include "Scene.h"
#include <cmath>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include "Pixel.h"
#include "tinyply.h"
#include "tinyxml2.h"
#define APPLY_FILTER_SINGLE_SAMPLE

inline float gaussian_filter(float x, float y, float sigma) {
  return exp(-(x * x + y * y) / (2 * sigma * sigma)) /
         (float)(2 * M_PI * sigma);
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
        Vector3 color = trace_ray(camera.calculate_ray_at(i + 0.5f, j + 0.5f),
                                  max_recursion_depth);
#ifdef APPLY_FILTER_SINGLE_SAMPLE
        for (int affected_j = j - 1; affected_j < j + 2; affected_j++) {
          if (affected_j < 0 || affected_j >= height) {
            continue;
          }
          for (int affected_i = i - 1; affected_i < i + 2; affected_i++) {
            if (affected_i < 0 || affected_i >= width) {
              continue;
            }
            float s_x = i - affected_i;
            float s_y = j - affected_j;
            result[affected_j * width + affected_i].add_color(
                color, gaussian_filter(s_x, s_y, 1.5f / 3.0f));
          }
        }
#else
        result[j * width + i].add_color(color, 1.0f);
#endif
      }
    }
  } else {
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<float> ms_distribution(0.0f, 1.0f);
    std::uniform_real_distribution<float> dof_distribution(-1.0f, 1.0f);
    std::uniform_real_distribution<float> time_distribution(0.0f, 1.0f);
    for (int j = starting_row; j < height; j += height_increase) {
      for (int i = 0; i < width; i++) {
        float aperture_size = camera.get_aperture_size();
        for (int x = 0; x < number_of_samples; x++) {
          for (int y = 0; y < number_of_samples; y++) {
            Vector3 color;
            float epsilon_x = ms_distribution(generator);
            float epsilon_y = ms_distribution(generator);
            float sample_x = (x + epsilon_x) / number_of_samples;
            float sample_y = (y + epsilon_y) / number_of_samples;
            if (aperture_size == 0.0f) {
              color = trace_ray(
                  camera.calculate_ray_at(i + sample_x, j + sample_y,
                                          time_distribution(generator)),
                  max_recursion_depth);
            } else {
              float dof_epsilon_x = dof_distribution(generator);
              float dof_epsilon_y = dof_distribution(generator);
              color =
                  trace_ray(camera.calculate_ray_at(
                                i + sample_x, j + sample_y, dof_epsilon_x,
                                dof_epsilon_y, time_distribution(generator)),
                            max_recursion_depth);
            }

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
                    color, gaussian_filter(s_x, s_y, 1.5f / 3.0f));
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
    if (ray.ray_type == r_primary) {
      if (background_texture) {
        return background_texture->get_color_at(ray.bg_u, ray.bg_v);
      } else {
        return background_color;
      }
    }
    return color;
  }
  const Shape* shape = hit_data.shape;
  const Vector3 intersection_point = ray.point_at(hit_data.t);
  const Material& material = materials[shape->get_material_id()];
  const int texture_id = shape->get_texture_id();
  const Texture* texture = (texture_id == -1) ? NULL : &(textures[texture_id]);
  const Vector3 w_0 = (ray.o - intersection_point).normalize();
  const Vector3& normal = hit_data.normal;

  Vector3 diffuse_color = material.diffuse;
  bool is_replace_all = false;
  if (texture) {
    is_replace_all = texture->get_decal_mode() == Texture::dm_replace_all;
    if (is_replace_all) {
      color = texture->get_color_at(hit_data.u, hit_data.v);
    } else {
      if (texture->is_perlin_noise()) {
        float perlin_value = hit_data.perlin_value;
        diffuse_color = texture->blend_color(
            Vector3(perlin_value, perlin_value, perlin_value), diffuse_color);
      } else {
        Vector3 texture_color = texture->get_color_at(hit_data.u, hit_data.v) /
                                texture->get_normalizer();
        diffuse_color = texture->blend_color(texture_color, diffuse_color);
      }
    }
  }
  if (!is_replace_all) {
    if (!ray.in_medium) {
      // ambient light
      // Should it be outside of if statement?
      color += material.ambient * ambient_light;

      // point lights
      for (const Light* light : lights) {
        // Shadow check
        bool has_diffuse = material.diffuse != zero_vector;
        bool has_specular = material.specular != zero_vector;
        if (!has_diffuse && !has_specular) {
          continue;
        }
        float light_distance;
        const Vector3 light_direction_vec =
            light->direction_and_distance(intersection_point, light_distance);
        const Vector3 w_i = light_direction_vec.normalize();
        Ray shadow_ray(intersection_point + (shadow_ray_epsilon * w_i), w_i,
                       r_shadow, ray.time);
        Hit_data shadow_hit_data;
        shadow_hit_data.t = std::numeric_limits<float>::infinity();
        shadow_hit_data.shape = NULL;
        bvh->intersect(shadow_ray, shadow_hit_data);
        if (shadow_hit_data.t < (light_distance - shadow_ray_epsilon) &&
            shadow_hit_data.t > 0.0f) {
          continue;
        }
        //
        Vector3 incoming_radiance =
            light->incoming_radiance(light_direction_vec);
        if (has_diffuse) {
          float diffuse_cos_theta = normal.dot(w_i);
          color += diffuse_color * incoming_radiance * diffuse_cos_theta;
        }
        if (has_specular) {
          float specular_cos_theta =
              fmax(0.0f, normal.dot((w_0 + w_i).normalize()));
          color += material.specular * incoming_radiance *
                   pow(specular_cos_theta, material.phong_exponent);
        }
      }
    }
    // Reflection
    if (material.mirror != zero_vector && current_recursion_depth > 0) {
      if (material.roughness == 0.0f) {
        const Vector3 w_r = ((2 * normal.dot(w_0) * normal) - w_0).normalize();
        Ray mirror_ray(intersection_point + (w_r * shadow_ray_epsilon), w_r,
                       r_reflection, ray.time);
        color += material.mirror *
                 trace_ray(mirror_ray, current_recursion_depth - 1);
      } else {
        const Vector3 w_r = ((2 * normal.dot(w_0) * normal) - w_0).normalize();
        Vector3 r_prime;
        if (w_r.x < w_r.y && w_r.x < w_r.z) {
          r_prime = Vector3(1.0f, w_r.y, w_r.z);
        } else if (w_r.y < w_r.z && w_r.y < w_r.x) {
          r_prime = Vector3(w_r.x, 1.0f, w_r.z);
        } else {
          r_prime = Vector3(w_r.x, w_r.y, 1.0f);
        }
        //(u,w_r,v) basis
        Vector3 u = r_prime.cross(w_r).normalize();
        Vector3 v = u.cross(w_r).normalize();
        std::mt19937 glossy_generator_u;
        glossy_generator_u.seed(
            std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_real_distribution<float> glossy_distribution_u(-0.5f,
                                                                    0.5f);
        std::mt19937 glossy_generator_v;
        glossy_generator_v.seed(
            std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_real_distribution<float> glossy_distribution_v(-0.5f,
                                                                    0.5f);
        float epsilon_u = glossy_distribution_u(glossy_generator_u);
        float epsilon_v = glossy_distribution_v(glossy_generator_v);
        const Vector3 w_r_prime =
            (w_r + material.roughness * (u * epsilon_u + v * epsilon_v))
                .normalize();
        Ray mirror_ray(intersection_point + (w_r_prime * shadow_ray_epsilon),
                       w_r_prime, r_reflection, ray.time);
        color += material.mirror *
                 trace_ray(mirror_ray, current_recursion_depth - 1);
      }
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
      if (d_n.dot(normal) < 0.0f) {
        refract_ray(d_n, normal, n, transmission_direction);
        cos_theta = (-d_n).dot(normal);
        k = Vector3(1.0f);
        entering_ray = true;
      } else {
        const Vector3& transparency = material.transparency;
        const float hit_data_t = hit_data.t;
        k.x = exp(log(transparency.x) * hit_data_t);
        k.y = exp(log(transparency.y) * hit_data_t);
        k.z = exp(log(transparency.z) * hit_data_t);
        entering_ray = false;
        if (refract_ray(d_n, -normal, 1.0f / n, transmission_direction)) {
          cos_theta = transmission_direction.dot(normal);
        } else {
          total_internal_reflection = true;
        }
      }

      if (total_internal_reflection) {
        Ray reflection_ray(intersection_point + (w_r * shadow_ray_epsilon), w_r,
                           r_reflection, ray.time);
        reflection_ray.in_medium = true;
        color += k * trace_ray(reflection_ray, current_recursion_depth - 1);
      } else {
        float r_0 = ((n - 1) * (n - 1)) / ((n + 1) * (n + 1));
        float r = r_0 + (1 - r_0) * pow(1.0f - cos_theta, 5);
        Ray reflection_ray(intersection_point + (w_r * shadow_ray_epsilon), w_r,
                           r_reflection, ray.time);
        reflection_ray.in_medium = !entering_ray;
        Ray transmission_ray(
            intersection_point + (transmission_direction * shadow_ray_epsilon),
            transmission_direction, r_refraction, ray.time);
        transmission_ray.in_medium = entering_ray;
        color +=
            k * (r * trace_ray(reflection_ray, current_recursion_depth - 1) +
                 (1 - r) *
                     trace_ray(transmission_ray, current_recursion_depth - 1));
      }
    }
  }
  return color;
}

Scene::Scene(const std::string& file_name) {
  const float degrees_to_radians = M_PI / 180.0f;
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
  debug("BackgroundColor is parsed");
  //

  // Get BackgroundTexture
  element = root->FirstChildElement("BackgroundTexture");
  if (element) {
    const std::string& bg_tex_name = element->GetText();
    background_texture = new Texture(bg_tex_name, "bilinear", "replace_all",
                                     "clamp", 255, 1.0f, false, 1.0f);
  } else {
    background_texture = nullptr;
  }
  debug("BackgroundTexture is parsed");
  //

  // Get ShadowRayEpsilon
  element = root->FirstChildElement("ShadowRayEpsilon");
  if (element) {
    stream << element->GetText() << std::endl;
  } else {
    stream << "0.001" << std::endl;
  }
  stream >> shadow_ray_epsilon;
  debug("ShadowRayEpsilon is parsed");
  //
  // Get IntersectionTestEpsilon
  /*element = root->FirstChildElement("IntersectionTestEpsilon");
  if (element) {
    stream << element->GetText() << std::endl;
    stream >> intersection_test_epsilon;
    debug("IntersectionTestEpsilon is parsed");
  }*/
  //
  // Get MaxRecursionDepth
  element = root->FirstChildElement("MaxRecursionDepth");
  if (element) {
    stream << element->GetText() << std::endl;
  } else {
    stream << "0" << std::endl;
  }
  stream >> max_recursion_depth;
  debug("MaxRecursionDepth is parsed");
  //

  // Get Cameras
  element = root->FirstChildElement("Cameras");
  element = element->FirstChildElement("Camera");
  while (element) {
    Vector3 position;
    Vector3 up;
    Vector3 gaze;
    float near_distance;
    float near_l, near_r, near_b, near_t;
    float focus_distance, aperture_size;
    int image_width, image_height;
    int number_of_samples;
    std::string image_name;

    auto child = element->FirstChildElement("Position");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("Up");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("NearDistance");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("FocusDistance");
    if (child) {
      stream << child->GetText() << std::endl;
    } else {
      stream << 0 << std::endl;
    }
    child = element->FirstChildElement("ApertureSize");
    if (child) {
      stream << child->GetText() << std::endl;
    } else {
      stream << 0.0f << std::endl;
    }
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
    stream >> focus_distance;
    stream >> aperture_size;
    stream >> image_width >> image_height;
    stream >> number_of_samples;
    number_of_samples = (int)sqrt(number_of_samples);
    if (number_of_samples <= 0) number_of_samples = 1;
    stream >> image_name;

    const char* camera_type = element->Attribute("type");
    if (camera_type && std::string(camera_type) == std::string("simple")) {
      child = element->FirstChildElement("GazePoint");
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
        float image_key, saturation_percentage, saturation, gamma;
        tonemap_child = element->FirstChildElement("TMOOptions");
        stream << tonemap_child->GetText() << std::endl;
        tonemap_child = element->FirstChildElement("Saturation");
        stream << tonemap_child->GetText() << std::endl;
        tonemap_child = element->FirstChildElement("Gamma");
        stream << tonemap_child->GetText() << std::endl;
        stream >> image_key >> saturation_percentage >> saturation >> gamma;
        tmo = new Photographic_tmo(image_key, saturation_percentage, saturation,
                                   gamma);
      }
    }
    if (tmo == nullptr) {
      tmo = new Null_tmo();
    }
    cameras.push_back(std::move(
        Camera(up, gaze, position, number_of_samples, image_name, near_l,
               near_r, near_b, near_t, near_distance, image_width, image_height,
               focus_distance, aperture_size, tmo)));
    element = element->NextSiblingElement("Camera");
  }
  stream.clear();
  debug("Cameras are parsed");
  // Cameras End

  // Get Lights
  element = root->FirstChildElement("Lights");
  auto child = element->FirstChildElement("AmbientLight");
  stream << child->GetText() << std::endl;
  stream >> ambient_light.x >> ambient_light.y >> ambient_light.z;
  element = element->FirstChildElement("PointLight");
  while (element) {
    Vector3 position;
    Vector3 intensity;
    child = element->FirstChildElement("Position");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("Intensity");
    stream << child->GetText() << std::endl;

    stream >> position.x >> position.y >> position.z;
    stream >> intensity.x >> intensity.y >> intensity.z;

    lights.push_back(new Point_light(position, intensity));
    element = element->NextSiblingElement("PointLight");
  }
  stream.clear();
  element = root->FirstChildElement("Lights");
  element = element->FirstChildElement("AreaLight");

  while (element) {
    Vector3 position;
    Vector3 intensity;
    Vector3 edge_vector_1, edge_vector_2;
    child = element->FirstChildElement("Position");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("Intensity");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("EdgeVector1");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("EdgeVector2");
    stream << child->GetText() << std::endl;
    stream >> position.x >> position.y >> position.z;
    stream >> intensity.x >> intensity.y >> intensity.z;
    stream >> edge_vector_1.x >> edge_vector_1.y >> edge_vector_1.z;
    stream >> edge_vector_2.x >> edge_vector_2.y >> edge_vector_2.z;
    lights.push_back(
        new Area_light(position, intensity, edge_vector_1, edge_vector_2));
    element = element->NextSiblingElement("AreaLight");
  }
  stream.clear();

  element = root->FirstChildElement("Lights");
  element = element->FirstChildElement("SpotLight");
  while (element) {
    Vector3 position;
    Vector3 intensity;
    Vector3 direction;
    float coverage_angle_in_radians, falloff_angle_in_radians;
    child = element->FirstChildElement("Position");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("Intensity");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("Direction");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("CoverageAngle");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("FalloffAngle");
    stream << child->GetText() << std::endl;
    stream >> position.x >> position.y >> position.z;
    stream >> intensity.x >> intensity.y >> intensity.z;
    stream >> direction.x >> direction.y >> direction.z;
    stream >> coverage_angle_in_radians >> falloff_angle_in_radians;
    coverage_angle_in_radians *= degrees_to_radians;
    falloff_angle_in_radians *= degrees_to_radians;
    lights.push_back(new Spot_light(position, intensity, direction,
                                    coverage_angle_in_radians,
                                    falloff_angle_in_radians));
    element = element->FirstChildElement("SpotLight");
  }
  stream.clear();

  element = root->FirstChildElement("Lights");
  element = element->FirstChildElement("DirectionalLight");
  while (element) {
    Vector3 direction;
    Vector3 radiance;
    float coverage_angle_in_radians, falloff_angle_in_radians;
    child = element->FirstChildElement("Direction");
    stream << child->GetText() << std::endl;
    child = element->FirstChildElement("Radiance");
    stream << child->GetText() << std::endl;
    stream >> direction.x >> direction.y >> direction.z;
    stream >> radiance.x >> radiance.y >> radiance.z;
    lights.push_back(new Directional_light(direction, radiance));
    element = element->FirstChildElement("DirectionalLight");
  }
  stream.clear();
  debug("Lights are parsed");
  // Lights End

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
    child = element->FirstChildElement("Roughness");
    if (child) {
      stream << child->GetText() << std::endl;
    } else {
      stream << "0" << std::endl;
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
    stream >> material.roughness;
    stream >> material.phong_exponent;
    stream >> material.transparency.x >> material.transparency.y >>
        material.transparency.z;
    stream >> material.refraction_index;
    if (element->BoolAttribute("degamma", false)) {
      material.ambient.x = pow(material.ambient.x, 0.45);
      material.ambient.y = pow(material.ambient.y, 0.45);
      material.ambient.z = pow(material.ambient.z, 0.45);
      material.diffuse.x = pow(material.diffuse.x, 0.45);
      material.diffuse.y = pow(material.diffuse.y, 0.45);
      material.diffuse.z = pow(material.diffuse.z, 0.45);
      material.specular.x = pow(material.specular.x, 0.45);
      material.specular.y = pow(material.specular.y, 0.45);
      material.specular.z = pow(material.specular.z, 0.45);
    }
    materials.push_back(material);
    element = element->NextSiblingElement("Material");
  }
  stream.clear();
  debug("Materials are parsed");
  // Materials end

  // Get Transformations
  element = root->FirstChildElement("Transformations");

  if (element) {
    // Get Scalings
    child = element->FirstChildElement("Scaling");
    while (child) {
      float x, y, z;
      stream << child->GetText() << std::endl;
      stream >> x >> y >> z;
      // TODO: maybe move?
      scaling_transformations.push_back(Scaling(x, y, z));
      child = child->NextSiblingElement("Scaling");
    }
    // Get Translations
    child = element->FirstChildElement("Translation");
    while (child) {
      float x, y, z;
      stream << child->GetText() << std::endl;
      stream >> x >> y >> z;
      // TODO: maybe move?
      translation_transformations.push_back(Translation(x, y, z));
      child = child->NextSiblingElement("Translation");
    }
    // Get Rotations
    child = element->FirstChildElement("Rotation");
    while (child) {
      float angle, x, y, z;
      stream << child->GetText() << std::endl;
      stream >> angle >> x >> y >> z;
      // TODO: maybe move?
      rotation_transformations.push_back(
          Rotation(angle * degrees_to_radians, x, y, z));
      child = child->NextSiblingElement("Rotation");
    }
  }
  stream.clear();
  debug("Transformations are parsed");
  // Transformations End

  // Get VertexData
  element = root->FirstChildElement("VertexData");
  if (element) {
    stream << element->GetText() << std::endl;
    Vector3 vertex;
    while (!(stream >> vertex.x).eof()) {
      stream >> vertex.y >> vertex.z;
      vertex_data.push_back(Vertex(vertex));
    }
  }
  stream.clear();
  debug("VertexData is parsed");
  // VertexData End

  // Get TexCoordData
  element = root->FirstChildElement("TexCoordData");
  if (element) {
    stream << element->GetText() << std::endl;
    Vector3 uv;
    while (!(stream >> uv.x).eof()) {
      stream >> uv.y;
      texture_coord_data.push_back(uv);
    }
    stream.clear();
    debug("TexCoordData is parsed");
  }
  // TexCoordData end

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
    const char* shading_mode = element->Attribute("shadingMode");
    Triangle_shading_mode triangle_shading_mode = tsm_flat;
    if (shading_mode && std::string(shading_mode) == std::string("smooth")) {
      triangle_shading_mode = tsm_smooth;
    }

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
    Vector3 velocity(0.0f);
    child = element->FirstChildElement("MotionBlur");
    if (child) {
      stream << child->GetText() << std::endl;
      stream >> velocity.x >> velocity.y >> velocity.z;
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
    std::vector<Shape*> triangles;

    child = element->FirstChildElement("Faces");
    int vertex_offset = vertex_data.size();
    int texture_offset = child->IntAttribute("textureOffset", 0);
    const char* ply_file = child->Attribute("plyFile");
    if (ply_file) {
      parse_ply_tinyply(std::string(ply_file), vertex_data, triangles,
                        vertex_offset, texture_offset, material_id, texture_id,
                        triangle_shading_mode);
    } else {
      vertex_offset = child->IntAttribute("vertexOffset", 0);
      stream << child->GetText() << std::endl;
      int v0_id, v1_id, v2_id;

      while (!(stream >> v0_id).eof()) {
        stream >> v1_id >> v2_id;
        v0_id--;
        v1_id--;
        v2_id--;
        Mesh_triangle* triangle = new Mesh_triangle(
            this, v0_id, v1_id, v2_id, vertex_offset, texture_offset,
            material_id, texture_id, triangle_shading_mode);
        float area = triangle->get_surface_area();
        const Vector3& surface_normal = triangle->normal;
        vertex_data[triangle->index_0 + triangle->offset].add_vertex_normal(
            surface_normal, area);
        vertex_data[triangle->index_1 + triangle->offset].add_vertex_normal(
            surface_normal, area);
        vertex_data[triangle->index_2 + triangle->offset].add_vertex_normal(
            surface_normal, area);
        triangles.push_back(triangle);
      }
    }
    stream.clear();
    meshes.push_back(
        new Mesh(material_id, texture_id, triangles,
                 Arbitrary_transformation(arbitrary_transformation), velocity));
    element = element->NextSiblingElement("Mesh");
  }
  stream.clear();
  debug("Meshes are parsed");

  // Create base mesh instances
  for (Mesh* mesh : meshes) {
    objects.push_back(new Mesh_instance(mesh->get_material_id(),
                                        mesh->texture_id, mesh,
                                        mesh->base_transform, mesh->velocity));
  }
  debug("Created base_mesh_instances");

  // Get MeshInstances
  element = root->FirstChildElement("Objects");
  element = element->FirstChildElement("MeshInstance");
  while (element) {
    Mesh* base_mesh = meshes[element->IntAttribute("baseMeshId") - 1];
    child = element->FirstChildElement("Material");
    int material_id;
    stream << child->GetText() << std::endl;
    stream >> material_id;
    material_id--;

    Matrix4x4 arbitrary_transformation =
        base_mesh->base_transform.get_transformation_matrix();
    const char* reset_transform = element->Attribute("resetTransform");
    if (reset_transform &&
        std::string(reset_transform) == std::string("true")) {
      arbitrary_transformation.make_identity();
    }

    // TODO: Until finding an elegant way, different textures for different mesh
    // instances are not supported. Since bump_map and perlin_noise calculations
    // are done in mesh_triangle of base_mesh Normal textures may work
    int texture_id = base_mesh->texture_id;
    /*child = element->FirstChildElement("Texture");
    if (child) {
      stream << child->GetText() << std::endl;
      stream >> texture_id;
      texture_id--;
    }
    stream.clear();*/

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
    Vector3 velocity(0.0f);
    child = element->FirstChildElement("MotionBlur");
    if (child) {
      stream << child->GetText() << std::endl;
      stream >> velocity.x >> velocity.y >> velocity.z;
    }
    stream.clear();
    objects.push_back(new Mesh_instance(
        material_id, texture_id, base_mesh,
        Arbitrary_transformation(arbitrary_transformation), velocity));
    element = element->NextSiblingElement("MeshInstance");
  }
  stream.clear();
  debug("MeshInstances are parsed");

  // Get Triangles
  element = root->FirstChildElement("Objects");
  element = element->FirstChildElement("Triangle");
  while (element) {
    int material_id;
    child = element->FirstChildElement("Material");
    stream << child->GetText() << std::endl;
    stream >> material_id;
    material_id--;
    int texture_id = -1;
    child = element->FirstChildElement("Texture");
    if (child) {
      stream << child->GetText() << std::endl;
      stream >> texture_id;
      texture_id--;
    }
    stream.clear();
    child = element->FirstChildElement("Indices");
    stream << child->GetText() << std::endl;
    int v0_id, v1_id, v2_id;
    stream >> v0_id >> v1_id >> v2_id;
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
    // No motion blur support for primitive triangles, yet
    objects.push_back(new Triangle(
        this, v0_id - 1, v1_id - 1, v2_id - 1, 0, material_id, texture_id,
        Arbitrary_transformation(arbitrary_transformation)));
    element = element->NextSiblingElement("Triangle");
  }
  stream.clear();
  debug("Triangles are parsed");

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
    Vector3 velocity(0.0f);
    child = element->FirstChildElement("MotionBlur");
    if (child) {
      stream << child->GetText() << std::endl;
      stream >> velocity.x >> velocity.y >> velocity.z;
    }
    stream.clear();
    objects.push_back(new Sphere(
        this, center_of_sphere, radius, material_id, texture_id,
        Arbitrary_transformation(arbitrary_transformation), velocity));
    element = element->NextSiblingElement("Sphere");
  }
  stream.clear();
  debug("Spheres are parsed");
  // Get Textures
  element = root->FirstChildElement("Textures");
  if (element) {
    element = element->FirstChildElement("Texture");
    while (element) {
      std::string image_name, interpolation_type, decal_mode, appearance;
      child = element->FirstChildElement("ImageName");
      image_name = child->GetText();
      child = element->FirstChildElement("Interpolation");
      if (child) {
        interpolation_type = child->GetText();
      } else {
        interpolation_type = "bilinear";
      }
      child = element->FirstChildElement("DecalMode");
      if (child) {
        decal_mode = child->GetText();
      } else {
        decal_mode = "blend_kd";
      }
      child = element->FirstChildElement("Appearance");
      if (child) {
        appearance = child->GetText();
      } else {
        appearance = "clamp";
      }
      child = element->FirstChildElement("Normalizer");
      if (child) {
        stream << child->GetText() << std::endl;
      } else {
        stream << "255.0" << std::endl;
      }
      child = element->FirstChildElement("ScalingFactor");
      if (child) {
        stream << child->GetText() << std::endl;
      } else {
        stream << "1.0" << std::endl;
      }
      float bumpmap_multiplier =
          element->FloatAttribute("bumpmapMultiplier", 1.0f);
      bool is_bump = element->BoolAttribute("bumpmap", false);
      float normalizer, scaling_factor;
      stream >> normalizer >> scaling_factor;
      textures.push_back(std::move(
          Texture(image_name, interpolation_type, decal_mode, appearance,
                  normalizer, scaling_factor, is_bump, bumpmap_multiplier)));

      element = element->NextSiblingElement("Texture");
    }
  }
  stream.clear();
  debug("Textures are parsed");
  bvh = BVH::create_bvh(objects);
  // Finalize surface normals
  for (Vertex& vertex : vertex_data) {
    if (vertex.has_vertex_normal()) {
      vertex.finalize_normal();
    }
  }
}

void Scene::parse_ply_tinyply(std::string filename,
                              std::vector<Vertex>& vertices,
                              std::vector<Shape*>& mesh_triangles,
                              int vertex_offset, int texture_offset,
                              int material_id, int texture_id,
                              Triangle_shading_mode tsm) const {
  try {
    // Read the file and create a std::istringstream suitable
    // for the lib -- tinyply does not perform any file i/o.
    std::ifstream ss(filename, std::ios::binary);

    if (ss.fail()) {
      throw std::runtime_error("failed to open " + filename);
    }
    tinyply::PlyFile file;
    file.parse_header(ss);
    std::shared_ptr<tinyply::PlyData> ply_vertices, ply_vertice_normals,
        ply_face_normals, ply_faces;
    try {
      ply_vertices =
          file.request_properties_from_element("vertex", {"x", "y", "z"});
    } catch (const std::exception& e) {
      std::cerr << "tinyply exception: " << e.what() << std::endl;
    }
    // vertice normals
    try {
      ply_vertice_normals =
          file.request_properties_from_element("vertex", {"nx", "ny", "nz"});
    } catch (const std::exception& e) {
      std::cerr << "tinyply exception: " << e.what() << std::endl;
    }
    // face normals
    try {
      ply_face_normals =
          file.request_properties_from_element("face", {"nx", "ny", "nz"});
    } catch (const std::exception& e) {
      std::cerr << "tinyply exception: " << e.what() << std::endl;
    }
    try {
      ply_faces =
          file.request_properties_from_element("face", {"vertex_index"});
    } catch (const std::exception& e) {
      try {
        ply_faces =
            file.request_properties_from_element("face", {"vertex_indices"});
      } catch (const std::exception& e) {
        std::cerr << "tinyply exception: " << e.what() << std::endl;
      }
      std::cerr << "tinyply exception: " << e.what() << std::endl;
    }
    file.read(ss);
    std::vector<Vector3> face_normals;

    // Process vertices
    if (ply_vertices) {
      // std::cout << "Parsing vertices" << std::endl;
      const size_t numVerticesBytes = ply_vertices->buffer.size_bytes();
      if (ply_vertices->t == tinyply::Type::FLOAT32) {
        std::vector<float> verts(ply_vertices->count * 3);
        std::memcpy(verts.data(), ply_vertices->buffer.get(), numVerticesBytes);
        size_t index = 0, c = ply_vertices->count;
        for (; index < c; index++) {
          vertices.push_back(Vector3(verts[3 * index], verts[3 * index + 1],
                                     verts[3 * index + 2]));
        }
      } else if (ply_vertices->t == tinyply::Type::FLOAT64) {
        // possible precision loss
        std::vector<double> verts(ply_vertices->count * 3);
        std::memcpy(verts.data(), ply_vertices->buffer.get(), numVerticesBytes);
        size_t index = 0, c = ply_vertices->count;
        for (; index < c; index++) {
          vertices.push_back(Vector3(verts[3 * index], verts[3 * index + 1],
                                     verts[3 * index + 2]));
        }
      } else {
        throw std::runtime_error("Unknown vertice type");
      }
    } else {
      throw std::runtime_error(filename + "contains no vertices");
    }
    // process vertice normals
    bool has_vertice_normals = false;
    if (ply_vertice_normals) {
      // std::cout << "Parsing normals" << std::endl;
      const size_t numNormalsBytes = ply_vertice_normals->buffer.size_bytes();
      if (ply_vertice_normals->t == tinyply::Type::FLOAT32) {
        std::vector<float> verts(ply_vertice_normals->count * 3);
        std::memcpy(verts.data(), ply_vertice_normals->buffer.get(),
                    numNormalsBytes);
        size_t index = 0, c = ply_vertice_normals->count;
        for (; index < c; index++) {
          Vector3 normal(verts[3 * index], verts[3 * index + 1],
                         verts[3 * index + 2]);
          vertices[index + vertex_offset].add_vertex_normal(normal, 1.0f);
          vertices[index + vertex_offset].finalize_normal();
        }
        has_vertice_normals = true;
      } else if (ply_vertice_normals->t == tinyply::Type::FLOAT64) {
        // possible precision loss
        std::vector<double> verts(ply_vertice_normals->count * 3);
        std::memcpy(verts.data(), ply_vertice_normals->buffer.get(),
                    numNormalsBytes);
        size_t index = 0, c = ply_vertice_normals->count;
        for (; index < c; index++) {
          Vector3 normal(verts[3 * index], verts[3 * index + 1],
                         verts[3 * index + 2]);
          vertices[index + vertex_offset].add_vertex_normal(normal, 1.0f);
          vertices[index + vertex_offset].finalize_normal();
        }
        has_vertice_normals = true;
      } else {
        throw std::runtime_error("Unknown normal type");
      }
    }
    // process face normals
    bool has_face_normals = false;
    if (ply_face_normals) {
      const size_t numFaceNormalsBytes = ply_face_normals->buffer.size_bytes();
      if (ply_face_normals->t == tinyply::Type::FLOAT32) {
        std::vector<float> vectors(ply_face_normals->count * 3);
        std::memcpy(vectors.data(), ply_face_normals->buffer.get(),
                    numFaceNormalsBytes);
        size_t index = 0, c = ply_face_normals->count;
        for (; index < c; index++) {
          face_normals.push_back(Vector3(vectors[3 * index],
                                         vectors[3 * index + 1],
                                         vectors[3 * index + 2])
                                     .normalize());
        }
        has_face_normals = true;
      } else if (ply_vertice_normals->t == tinyply::Type::FLOAT64) {
        // Possible precision loss
        std::vector<double> vectors(ply_face_normals->count * 3);
        std::memcpy(vectors.data(), ply_face_normals->buffer.get(),
                    numFaceNormalsBytes);
        size_t index = 0, c = ply_face_normals->count;
        for (; index < c; index++) {
          face_normals.push_back(Vector3(vectors[3 * index],
                                         vectors[3 * index + 1],
                                         vectors[3 * index + 2])
                                     .normalize());
        }
        has_face_normals = true;
      } else {
        throw std::runtime_error("Unknown normal type");
      }
    }
    // process faces
    if (ply_faces) {
      // std::cout << "Parsing faces" << std::endl;
      const size_t numFacesBytes = ply_faces->buffer.size_bytes();
      size_t index = 0, c = ply_faces->count;
      if (ply_faces->t == tinyply::Type::INT32) {
        bool is_triangle = c * 12 == numFacesBytes;
        bool is_quad = c * 16 == numFacesBytes;
        if (sizeof(int) != 4) {
          throw std::runtime_error("sizeof(int)!=4");
        }
        std::vector<int> verts(numFacesBytes / 4);
        std::memcpy(verts.data(), ply_faces->buffer.get(), numFacesBytes);
        if (is_triangle) {
          for (; index < c; index++) {
            int index_0 = verts[index * 3];
            int index_1 = verts[index * 3];
            int index_2 = verts[index * 3];
            mesh_triangles.push_back(new Mesh_triangle(
                this, index_0, index_1, index_2, vertex_offset, texture_offset,
                material_id, texture_id, tsm));
            Mesh_triangle* triangle =
                (Mesh_triangle*)*(mesh_triangles.rbegin());
            if (has_face_normals) {
              triangle->normal = face_normals[index];
            }
            if (!has_vertice_normals) {
              float area = triangle->get_surface_area();
              const Vector3& surface_normal = triangle->normal;
              vertices[index_0 + vertex_offset].add_vertex_normal(
                  surface_normal, area);
              vertices[index_1 + vertex_offset].add_vertex_normal(
                  surface_normal, area);
              vertices[index_2 + vertex_offset].add_vertex_normal(
                  surface_normal, area);
            }
          }

        } else if (is_quad) {
          for (; index < c; index++) {
            int index_0 = verts[index * 4];
            int index_1 = verts[index * 4 + 1];
            int index_2 = verts[index * 4 + 2];
            int index_3 = verts[index * 4 + 3];

            mesh_triangles.push_back(new Mesh_triangle(
                this, index_0, index_1, index_3, vertex_offset, texture_offset,
                material_id, texture_id, tsm));
            mesh_triangles.push_back(new Mesh_triangle(
                this, index_2, index_3, index_1, vertex_offset, texture_offset,
                material_id, texture_id, tsm));
            Mesh_triangle* triangle1 =
                (Mesh_triangle*)*(mesh_triangles.rbegin()++);
            Mesh_triangle* triangle2 =
                (Mesh_triangle*)*(mesh_triangles.rbegin());
            if (has_face_normals) {
              triangle1->normal = face_normals[index];
              triangle2->normal = face_normals[index];
            }
            if (!has_vertice_normals) {
              {
                float area = triangle1->get_surface_area();
                const Vector3& surface_normal = triangle1->normal;
                vertices[index_0 + vertex_offset].add_vertex_normal(
                    surface_normal, area);
                vertices[index_1 + vertex_offset].add_vertex_normal(
                    surface_normal, area);
                vertices[index_2 + vertex_offset].add_vertex_normal(
                    surface_normal, area);
              }
              {
                float area = triangle2->get_surface_area();
                const Vector3& surface_normal = triangle2->normal;
                vertices[index_2 + vertex_offset].add_vertex_normal(
                    surface_normal, area);
                vertices[index_3 + vertex_offset].add_vertex_normal(
                    surface_normal, area);
                vertices[index_1 + vertex_offset].add_vertex_normal(
                    surface_normal, area);
              }
            }
          }
        } else {
          throw std::runtime_error("This parser doesn't support hybrid files");
        }
      } else {
        throw std::runtime_error("check vertex_index type type");
      }
    } else {
      throw std::runtime_error(filename + "contains no faces");
    }
  } catch (const std::exception& e) {
    std::cerr << "Caught tinyply exception: " << e.what() << std::endl;
    exit(-1);
  }
}
Scene::~Scene() {
  delete bvh;
  if (background_texture) {
    delete background_texture;
  }
  size_t size = meshes.size();
  for (size_t i = 0; i < size; i++) {
    delete meshes[i];
  }
  size = lights.size();
  for (size_t i = 0; i < size; i++) {
    delete lights[i];
  }
}
