#include "Scene.h"
#include <cmath>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include "Bounding_volume_hierarchy.h"
#include "Mesh.h"
#include "Point_light.h"
#include "Sphere.h"
#include "tinyxml2.h"
//#define GAUSSIAN_FILTER
#define ALPHA 0.7f
void Scene::reset_hash_grid() {
  hash_grid.clear();
  for (int i = 0; i < hit_points.size(); i++) {
    delete hit_points[i];
  }
  hit_points.clear();
}
void Scene::build_hash_grid(const int width, const int height) {
  hit_point_bbox = Bounding_box();
  for (int i = 0; i < hit_points.size(); i++) {
    Hit_point* hit_point = hit_points[i];
    hit_point_bbox.fit(hit_point->position);
  }
  Vector3 bbox_size = hit_point_bbox.delta;
  float initial_radius = ((bbox_size.x + bbox_size.y + bbox_size.z) / 3.0f) /
                         ((width + height) / 2.0f) * 2.0f * 4.0f;
  hit_point_bbox = Bounding_box();
  num_hash = hit_points.size();
  for (int i = 0; i < hit_points.size(); i++) {
    Hit_point* hit_point = hit_points[i];
    hit_point->radius_squared = initial_radius * initial_radius;
    hit_point->n = 0;
    hit_point->flux = Vector3(0.0f);
    hit_point_bbox.fit(hit_point->position - initial_radius);
    hit_point_bbox.fit(hit_point->position + initial_radius);
  }
  hash_scale = 1.0 / (initial_radius * 2.0);

  hash_grid.resize(num_hash);
  for (int i = 0; i < hit_points.size(); i++) {
    Hit_point* hit_point = hit_points[i];
    Vector3 BMin =
        ((hit_point->position - initial_radius) - hit_point_bbox.min_corner) *
        hash_scale;
    Vector3 BMax =
        ((hit_point->position + initial_radius) - hit_point_bbox.min_corner) *
        hash_scale;
    for (int iz = std::abs(int(BMin.z)); iz <= std::abs(int(BMax.z)); iz++) {
      for (int iy = std::abs(int(BMin.y)); iy <= std::abs(int(BMax.y)); iy++) {
        for (int ix = std::abs(int(BMin.x)); ix <= std::abs(int(BMax.x));
             ix++) {
          int hv = hash(ix, iy, iz);
          hash_grid[hv].push_back(hit_point);
        }
      }
    }
  }
}
void Scene::photon_trace(const Ray& ray, int depth, const Vector3& flux,
                         const Vector3& attenuation, int photon_id) {
  Intersection intersection;
  if (!bvh->intersect(ray, intersection, true)) {
    return;
  }
  int depth3 = (depth + 1) * 3;
  const Shape* shape = intersection.shape;
  Vector3 intersection_point = ray.point_at(intersection.t);
  const Vector3& normal = intersection.normal;
  Vector3 nl = normal.dot(ray.d) < 0 ? normal : normal * -1;
  const Material& material = materials[shape->get_material_id()];
  if (material.material_type == mt_diffuse) {
    // Use Quasi-Monte Carlo to sample the next direction
    float r1 = 2. * M_PI * Light::hal(depth3 - 1, photon_id);
    float r2 = Light::hal(depth3 + 0, photon_id);
    float r2s = std::sqrt(r2);

    Vector3 hh = (intersection_point - hit_point_bbox.min_corner) * hash_scale;
    int ix = abs(int(hh.x));
    int iy = abs(int(hh.y));
    int iz = abs(int(hh.z));
    {
      std::vector<Hit_point*>& hpoints = hash_grid[hash(ix, iy, iz)];
      for (int i = 0; i < hpoints.size(); i++) {
        Hit_point* hit_point = hpoints[i];

        Vector3 v = hit_point->position - intersection_point;
        hit_point->mutex.lock();
        if ((hit_point->normal.dot(normal) > 1e-3f) &&
            (v.dot(v) <= hit_point->radius_squared)) {
          // Unlike N in the paper, hit_point->n stores "N / ALPHA" to make it
          // an integer value
          float radius_reduction =
              (hit_point->n * ALPHA + ALPHA) / (hit_point->n * ALPHA + 1.0);
          hit_point->radius_squared =
              hit_point->radius_squared * radius_reduction;
          hit_point->n++;
          Vector3 color(0.0f);
          Vector3 w_i = -ray.d.normalize();
          if (hit_point->material.brdf_id == -1) {
            // all things except w_i should be used from hit_point
            float cos_theta_i = hit_point->normal.dot(w_i);
            if (cos_theta_i > 1.0f || cos_theta_i <= 0.0f) {
              color = 0.0f;
            } else {
              float specular_cos_theta = std::max(
                  0.0f,
                  hit_point->normal.dot((hit_point->w_o + w_i).normalize()));
              color = (hit_point->material.diffuse +
                       hit_point->material.specular *
                           std::pow(specular_cos_theta,
                                    hit_point->material.phong_exponent) /
                           cos_theta_i) *
                      hit_point->attenuation;
            }

          } else {
            // parse brdfs, use brdfs with hit_point's diffuse, specular
          }
          hit_point->flux = (hit_point->flux + color * flux * (1.0f / M_PI)) *
                            radius_reduction;
        }
        hit_point->mutex.unlock();
      }
    }
    Vector3 w = nl;
    Vector3 u = ((std::fabs(w.x) > 0.1f) ? Vector3(0.0f, 1.0f, 0.0f)
                                         : Vector3(1.0f, 0.0f, 0.0f))
                    .cross(w)

                    .normalize();
    Vector3 v = w.cross(u);
    Vector3 d = (u * std::cos(r1) * r2s + v * std::sin(r1) * r2s +
                 w * std::sqrt(1.0f - r2))
                    .normalize();
    Vector3 base_color(0.0f);
    Vector3 w_i = -ray.d.normalize();
    const Vector3& w_o = d;

    if (material.brdf_id == -1) {
      float cos_theta_i = nl.dot(w_i);
      if (cos_theta_i > 1.0f || cos_theta_i <= 0.0f) {
        base_color = 0.0f;
      } else {
        float specular_cos_theta =
            std::max(0.0f, nl.dot((w_o + w_i).normalize()));
        base_color = (material.diffuse + material.specular *
                                             std::pow(specular_cos_theta,
                                                      material.phong_exponent) /
                                             cos_theta_i) *
                     attenuation;
      }
    } else {
      // parse brdfs, use brdfs with intersection's diffuse, specular
    }
    float p = base_color.x > base_color.y && base_color.x > base_color.z
                  ? base_color.x
                  : base_color.y > base_color.z ? base_color.y : base_color.z;
    // Be sure that hal cannot return negative
    if (depth < max_recursion_depth && Light::hal(depth3 + 1, photon_id) < p) {
      // Ray_type is not important for ppm, not checking it
      photon_trace(
          Ray(intersection_point + (d * shadow_ray_epsilon), d, r_reflection),
          depth + 1, base_color * (flux) * (1.0f / p), attenuation, photon_id);
    }
  } else if (depth >= max_recursion_depth) {
    return;
  } else if (material.material_type == mt_mirror) {
    const Vector3 w_o = (ray.o - intersection_point).normalize();
    const Vector3 w_r = ((2.0f * normal.dot(w_o) * normal) - w_o).normalize();
    Ray mirror_ray(intersection_point + (w_r * shadow_ray_epsilon), w_r,
                   r_reflection);
    photon_trace(mirror_ray, depth + 1, material.mirror * flux,
                 material.mirror * attenuation, photon_id);

  } else if (material.material_type == mt_refractive) {
    // TODO: this calculation does not calculate attenuation with
    // exp(log(transparency) * hit_data_t).
    // uses only transparency while refracting ray, maybe we can convert it to
    // our formula, AKCAY?
    const Vector3 nl = normal.dot(ray.d) < 0.0f ? normal : normal * -1;
    const Vector3 w_o = (ray.o - intersection_point).normalize();
    const Vector3 w_r = ((2.0f * normal.dot(w_o) * normal) - w_o).normalize();
    Ray reflection_ray(intersection_point + (w_r * shadow_ray_epsilon), w_r,
                       r_reflection);
    bool into = (normal.dot(nl) > 0.0f);
    constexpr float air_index = 1.0f;
    float nnt = into ? air_index / material.refraction_index
                     : material.refraction_index / air_index;
    float ddn = ray.d.dot(nl);
    float cos2t = 1 - nnt * nnt * (1 - ddn * ddn);
    if (cos2t < 0.0f) {
      photon_trace(reflection_ray, depth + 1, flux, attenuation, photon_id);
      return;
    }
    Vector3 refraction_direction =
        (ray.d * nnt -
         normal * ((into ? 1 : -1) * (ddn * nnt + std::sqrt(cos2t))))
            .normalize();
    float a = material.refraction_index - air_index;
    float b = material.refraction_index + air_index;
    float R0 = a * a / (b * b);

    float cosinealpha = into ? -ddn : refraction_direction.dot(normal);
    float c = 1 - cosinealpha;
    float fresnel = R0 + (1 - R0) * c * c * c * c * c;
    float P = fresnel;
    Ray refraction_ray(
        intersection_point + (refraction_direction * shadow_ray_epsilon),
        refraction_direction, r_refraction);
    Vector3 attenuated_color = material.transparency * attenuation;
    if (Light::hal(depth3 - 1, photon_id) < P) {
      photon_trace(reflection_ray, depth + 1, flux, attenuated_color,
                   photon_id);
    } else {
      photon_trace(refraction_ray, depth + 1, flux, attenuated_color,
                   photon_id);
    }
  }
}
void Scene::eye_trace_lines(int index, int starting_row, int height_increase) {
  const Camera& camera = cameras[index];
  const Image_plane& image_plane = camera.get_image_plane();
  const int width = image_plane.width;
  const int height = image_plane.height;
  int number_of_samples = camera.get_number_of_samples();
  if (number_of_samples == 1) {
    for (int j = starting_row; j < height; j += height_increase) {
      for (int i = 0; i < width; i++) {
        Ray primary_ray = camera.calculate_ray_at(i + 0.5f, j + 0.5f);
        eye_trace(primary_ray, 0, 1.0f, j * width + i, 1.0f);
      }
    }
  } else {
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<float> ms_distribution(0.0f, 1.0f);
    // Maximum 4 samples for one pixel
    number_of_samples = std::min(2, number_of_samples);
    for (int j = starting_row; j < height; j += height_increase) {
      for (int i = 0; i < width; i++) {
        for (int x = 0; x < number_of_samples; x++) {
          for (int y = 0; y < number_of_samples; y++) {
            float epsilon_x = ms_distribution(generator);
            float epsilon_y = ms_distribution(generator);
            float sample_x = (x + epsilon_x) / number_of_samples;
            float sample_y = (y + epsilon_y) / number_of_samples;
            Ray primary_ray =
                camera.calculate_ray_at(i + sample_x, j + sample_y);
            eye_trace(primary_ray, 0, 1.0f, j * width + i, 1.0f);
          }
        }
      }
    }
  }
}
void Scene::eye_trace(const Ray& ray, int depth, const Vector3& attenuation,
                      unsigned int pixel_index, float pixel_weight) {
  Intersection intersection;
  if (!bvh->intersect(ray, intersection, true)) {
    return;
  }
  const Shape* shape = intersection.shape;
  Vector3 intersection_point = ray.point_at(intersection.t);
  const Vector3& normal = intersection.normal;
  const Material& material = materials[shape->get_material_id()];
  if (material.material_type == mt_diffuse) {
    Hit_point* hit_point = new Hit_point();
    hit_point->material = material;
    hit_point->attenuation = attenuation;
    hit_point->w_o = (ray.o - intersection_point).normalize();
    hit_point->normal = normal;
    hit_point->position = intersection_point;
    hit_point->pixel = pixel_index;
    hit_point->pixel_weight = pixel_weight;
    add_hit_point(hit_point);
  } else if (depth >= max_recursion_depth) {
    return;
  } else if (material.material_type == mt_mirror) {
    const Vector3 w_o = (ray.o - intersection_point).normalize();
    const Vector3 w_r = ((2.0f * normal.dot(w_o) * normal) - w_o).normalize();
    Ray mirror_ray(intersection_point + (w_r * shadow_ray_epsilon), w_r,
                   r_reflection);
    eye_trace(mirror_ray, depth + 1, material.mirror * attenuation, pixel_index,
              pixel_weight);
  } else if (material.material_type == mt_refractive) {
    // TODO: this calculation does not calculate attenuation with
    // exp(log(transparency) * hit_data_t).
    // uses only transparency while refracting ray, maybe we can convert it to
    // our formula, AKCAY?
    const Vector3 nl = normal.dot(ray.d) < 0.0f ? normal : normal * -1;
    const Vector3 w_o = (ray.o - intersection_point).normalize();
    const Vector3 w_r = ((2.0f * normal.dot(w_o) * normal) - w_o).normalize();
    Ray reflection_ray(intersection_point + (w_r * shadow_ray_epsilon), w_r,
                       r_reflection);
    bool into = (normal.dot(nl) > 0.0f);
    constexpr float air_index = 1.0f;
    float nnt = into ? air_index / material.refraction_index
                     : material.refraction_index / air_index;
    float ddn = ray.d.dot(nl);
    float cos2t = 1 - nnt * nnt * (1 - ddn * ddn);
    if (cos2t < 0.0f) {
      eye_trace(reflection_ray, depth + 1, attenuation, pixel_index,
                pixel_weight);
      return;
    }
    Vector3 refraction_direction =
        (ray.d * nnt -
         normal * ((into ? 1 : -1) * (ddn * nnt + std::sqrt(cos2t))))
            .normalize();
    float a = material.refraction_index - air_index;
    float b = material.refraction_index + air_index;
    float R0 = a * a / (b * b);

    float cosinealpha = into ? -ddn : refraction_direction.dot(normal);
    float c = 1 - cosinealpha;
    float fresnel = R0 + (1 - R0) * c * c * c * c * c;
    Ray refraction_ray(
        intersection_point + (refraction_direction * shadow_ray_epsilon),
        refraction_direction, r_refraction);
    Vector3 attenuated_color = material.transparency * attenuation;
    eye_trace(reflection_ray, depth + 1, attenuated_color * fresnel,
              pixel_index, pixel_weight);
    eye_trace(refraction_ray, depth + 1, attenuated_color * (1.0 - fresnel),
              pixel_index, pixel_weight);
  }
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
  // Get Lights
  element = root->FirstChildElement("Lights");
  if (element) {
    Point_light::load_point_lights_from_xml(element, lights);
  }
  //
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
