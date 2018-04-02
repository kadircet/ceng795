#include "Scene.h"
#include <cmath>
#include <random>
#include <sstream>
#include <string>
#include "Pixel.h"
#include "tinyxml2.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

inline float gaussian_filter(float x, float y, float sigma) {
  return exp(-(x * x + y * y) / (2 * sigma * sigma)) / (float) (2 * M_PI * sigma);
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
        std::default_random_engine ms_generator;
        ms_generator.seed(
            std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_real_distribution<float> ms_distribution(0.0, 1);
        float aperture_size = camera.get_aperture_size();
        std::default_random_engine dof_generator;
        dof_generator.seed(
          std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_real_distribution<float> dof_distribution(-1.0f, 1.0f);
        for (int x = 0; x < number_of_samples; x++) {
          for (int y = 0; y < number_of_samples; y++) {
            Vector3 color;
            float epsilon_x = ms_distribution(ms_generator);
            float epsilon_y = ms_distribution(ms_generator);
            float sample_x = (x + epsilon_x) / number_of_samples;
            float sample_y = (y + epsilon_y) / number_of_samples;
            if (aperture_size == 0.0f) {
              // since calculate_ray_at adds 0.5 to the pixel number, subtracting
              // 0.5
              color = trace_ray(
                camera.calculate_ray_at(i + sample_x - 0.5, j + sample_y - 0.5),
                max_recursion_depth);
            }
            else {
              float dof_epsilon_x = dof_distribution(dof_generator);
              float dof_epsilon_y = dof_distribution(dof_generator);
              // since calculate_ray_at adds 0.5 to the pixel number, subtracting
              // 0.5
              color = trace_ray(
                camera.calculate_ray_at(i + sample_x - 0.5, j + sample_y - 0.5, dof_epsilon_x, dof_epsilon_y),
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
      if(material.diffuse != zero_vector) {
        float diffuse_cos_theta = normal.dot(w_i);
        color += material.diffuse * point_light.intensity * diffuse_cos_theta /
                 light_distance_squared;
      }
      if(material.specular != zero_vector) {
        float specular_cos_theta =
          fmax(0.0f, normal.dot((w_0 + w_i).normalize()));
        color += material.specular * point_light.intensity *
          pow(specular_cos_theta, material.phong_exponent) /
          light_distance_squared;
      }
    }
  }
  // Reflection
  if (material.mirror != zero_vector && current_recursion_depth > 0) {
    if (material.roughness == 0.0f) {
      const Vector3 w_r = ((2 * normal.dot(w_0) * normal) - w_0).normalize();
      Ray mirror_ray(intersection_point + (w_r * shadow_ray_epsilon), w_r);
      color +=
        material.mirror * trace_ray(mirror_ray, current_recursion_depth - 1);
    } else {
      /*const Vector3 w_r = ((2 * normal.dot(w_0) * normal) - w_0).normalize();
      Vector3 t;
      float x = std::fabs(w_r.x);
      float y = std::fabs(w_r.y);
      float z = std::fabs(w_r.z);
      if (x < y && x < z) {
        t = Vector3(1.0f, w_r.y, w_r.z);
      }
      else if (y < z && y < x) {
        t = Vector3(w_r.x, 1.0f, w_r.z);
      }
      else {
        t = Vector3(w_r.x, w_r.y, 1.0f);
      }
      //(u,w_r,v) basis
      Vector3 u = t.cross(w_r).normalize();
      Vector3 v = u.cross(w_r).normalize();
      std::default_random_engine glossy_generator;
      glossy_generator.seed(
        std::chrono::system_clock::now().time_since_epoch().count());
      std::uniform_real_distribution<float> glossy_distribution(-1.0f, 1.0f);
      float epsilon_u = glossy_distribution(glossy_generator);
      float epsilon_v = glossy_distribution(glossy_generator);
      const Vector3 w_r_prime = (w_r + material.roughness*(u*epsilon_u + v * epsilon_v)).normalize();
      Ray mirror_ray(intersection_point + (w_r_prime * shadow_ray_epsilon), w_r_prime);
      color +=
        material.mirror * trace_ray(mirror_ray, current_recursion_depth - 1);*/
      const Vector3 w = ((2 * normal.dot(w_0) * normal) - w_0).normalize();
      Vector3 t;
      float x = std::fabs(w.x);
      float y = std::fabs(w.y);
      float z = std::fabs(w.z);
      if (x < y && x < z) {
        t = Vector3(1.0f, w.y, w.z);
      }
      else if (y < z && y < x) {
        t = Vector3(w.x, 1.0f, w.z);
      }
      else {
        t = Vector3(w.x, w.y, 1.0f);
      }
      const Vector3 u = t.cross(w).normalize();
      const Vector3 v = w.cross(u);
      std::default_random_engine glossy_generator;
      glossy_generator.seed(
        std::chrono::system_clock::now().time_since_epoch().count());
      std::uniform_real_distribution<float> glossy_distribution(0.0f, 1.0f);
      float epsilon_u = glossy_distribution(glossy_generator);
      float epsilon_v = glossy_distribution(glossy_generator);
      const Vector3 r_prime = ((epsilon_u - 0.5)*material.roughness * u + (epsilon_v - 0.5)*material.roughness * v + w).normalize();
      Ray mirror_ray(intersection_point + (r_prime * shadow_ray_epsilon), r_prime);
      color +=
        material.mirror * trace_ray(mirror_ray, current_recursion_depth - 1);
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
      Ray reflection_ray(intersection_point + (w_r * shadow_ray_epsilon), w_r);
      reflection_ray.in_medium = true;
      color += k*trace_ray(reflection_ray, current_recursion_depth - 1);
    } else {
      float r_0 = ((n - 1) * (n - 1)) / ((n + 1) * (n + 1));
      float r = r_0 + (1 - r_0) * pow(1.0f - cos_theta, 5);
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
  debug("BackgroundColor is parsed");

  // Get ShadowRayEpsilon
  element = root->FirstChildElement("ShadowRayEpsilon");
  if (element) {
    stream << element->GetText() << std::endl;
  } else {
    stream << "0.001" << std::endl;
  }
  stream >> shadow_ray_epsilon;
  debug("ShadowRayEpsilon is parsed");

  // Get MaxRecursionDepth
  element = root->FirstChildElement("MaxRecursionDepth");
  if (element) {
    stream << element->GetText() << std::endl;
  } else {
    stream << "0" << std::endl;
  }
  stream >> max_recursion_depth;
  debug("MaxRecursionDepth is parsed");

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
    child = element->FirstChildElement("FocusDistance");
    if (child) {
      stream << child->GetText() << std::endl;
    }
    else {
      stream << 0 << std::endl;
    }
    child = element->FirstChildElement("ApertureSize");
    if (child) {
      stream << child->GetText() << std::endl;
    }
    else {
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
    Vector3 position, up, gaze;
    float near_distance;
    float near_l, near_r, near_b, near_t;
    float focus_distance, aperture_size;
    int image_width, image_height;
    int number_of_samples;
    std::string image_name;
    stream >> position.x >> position.y >> position.z;
    stream >> gaze.x >> gaze.y >> gaze.z;
    stream >> up.x >> up.y >> up.z;
    stream >> near_l >> near_r >> near_b >> near_t;
    stream >> near_distance;
    stream >> focus_distance;
    stream >> aperture_size;
    stream >> image_width >> image_height;
    stream >> number_of_samples;
    number_of_samples =(int) sqrt(number_of_samples);
    if (number_of_samples <= 0) number_of_samples = 1;
    stream >> image_name;
    Camera camera(up, gaze, position, number_of_samples, image_name, near_l,
                  near_r, near_b, near_t, near_distance, image_width,
                  image_height,focus_distance,aperture_size);
    cameras.push_back(camera);
    element = element->NextSiblingElement("Camera");
  }
  stream.clear();
  debug("Cameras are parsed");

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
  stream.clear();
  debug("Lights are parsed");

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

    materials.push_back(material);
    element = element->NextSiblingElement("Material");
  }
  stream.clear();
  debug("Materials are parsed");

  // Get Transformations
  element = root->FirstChildElement("Transformations");
  const float degree_to_pi = M_PI / 180.0f;
  if (element) {
	  // Get Scalings
	  child = element->FirstChildElement("Scaling");
	  while (child) {
		  float x, y, z;
		  stream << child->GetText() << std::endl;
		  stream >> x >> y >> z;
		  //TODO: maybe move?
		  scaling_transformations.push_back(Scaling(x, y, z));
		  child = child->NextSiblingElement("Scaling");
	  }
	  // Get Translations
	  child = element->FirstChildElement("Translation");
	  while (child) {
		  float x, y, z;
		  stream << child->GetText() << std::endl;
		  stream >> x >> y >> z;
		  //TODO: maybe move?
		  translation_transformations.push_back(Translation(x, y, z));
		  child = child->NextSiblingElement("Translation");
	  }
	  // Get Rotations
	  child = element->FirstChildElement("Rotation");
	  while (child) {
		  float angle, x, y, z;
		  stream << child->GetText() << std::endl;
		  stream >> angle >> x >> y >> z;
		  //TODO: maybe move?
		  rotation_transformations.push_back(Rotation(angle*degree_to_pi, x, y, z));
		  child = child->NextSiblingElement("Rotation");
	  }
  }
  stream.clear();
  debug("Transformations are parsed");

  // Get VertexData
  element = root->FirstChildElement("VertexData");
  stream << element->GetText() << std::endl;
  Vector3 vertex;
  while (!(stream >> vertex.x).eof()) {
    stream >> vertex.y >> vertex.z;
    vertex_data.push_back(Vertex(vertex));
  }
  stream.clear();
  debug("VertexData is parsed");


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
    if(shading_mode && std::string(shading_mode) == std::string("smooth"))
    {
      triangle_shading_mode=tsm_smooth;
    }
    child = element->FirstChildElement("Faces");
    int vertex_offset = child->IntAttribute("vertexOffset",0);
    stream << child->GetText() << std::endl;
    int v0_id, v1_id, v2_id;
    std::vector<Shape*> triangles;
    while (!(stream >> v0_id).eof()) {
      stream >> v1_id >> v2_id;
      triangles.push_back(
          new Mesh_triangle(this, v0_id - 1, v1_id - 1, v2_id - 1, vertex_offset, material_id,triangle_shading_mode));
    }
    stream.clear();

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
				  arbitrary_transformation = scaling_transformations[index].get_transformation_matrix() * arbitrary_transformation;
				  break;
			  case 't':
				  arbitrary_transformation = translation_transformations[index].get_transformation_matrix() * arbitrary_transformation;
				  break;
			  case 'r':
				  arbitrary_transformation = rotation_transformations[index].get_transformation_matrix() * arbitrary_transformation;
				  break;
			  }
		  }
		  stream.clear();
	  }
    meshes.push_back(new Mesh(material_id, -1, triangles, Arbitrary_transformation(arbitrary_transformation)));
    
    //Calculate vertex normals
    for(Shape* shape: triangles) {
      Mesh_triangle* triangle = (Mesh_triangle*) shape;
      float area = triangle->get_surface_area();
      const Vector3& surface_normal = triangle->normal;
      vertex_data[triangle->index_0+triangle->offset].add_vertex_normal(surface_normal, area);
      vertex_data[triangle->index_1+triangle->offset].add_vertex_normal(surface_normal, area);
      vertex_data[triangle->index_2+triangle->offset].add_vertex_normal(surface_normal, area);
    }
    element = element->NextSiblingElement("Mesh");
  }
  stream.clear();
  debug("Meshes are parsed");
  
  //Create base mesh instances
  for (Mesh* mesh : meshes)
  {
	  objects.push_back(new Mesh_instance(mesh->get_material_id(), mesh->texture_id, mesh, mesh->base_transform));
  }

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
    Matrix4x4 arbitrary_transformation(true);
    //TODO: Check hard reset parameter
    if (true) {
      arbitrary_transformation = base_mesh->base_transform.get_transformation_matrix();
    }
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
          arbitrary_transformation = scaling_transformations[index].get_transformation_matrix() * arbitrary_transformation;
          break;
        case 't':
          arbitrary_transformation = translation_transformations[index].get_transformation_matrix() * arbitrary_transformation;
          break;
        case 'r':
          arbitrary_transformation = rotation_transformations[index].get_transformation_matrix() * arbitrary_transformation;
          break;
        }
      }
      stream.clear();
    }
    objects.push_back(new Mesh_instance(material_id, -1, base_mesh, Arbitrary_transformation(arbitrary_transformation)));
    element = element->NextSiblingElement("MeshInstance");
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
				  arbitrary_transformation = scaling_transformations[index].get_transformation_matrix() * arbitrary_transformation;
				  break;
			  case 't':
				  arbitrary_transformation = translation_transformations[index].get_transformation_matrix() * arbitrary_transformation;
				  break;
			  case 'r':
				  arbitrary_transformation = rotation_transformations[index].get_transformation_matrix() * arbitrary_transformation;
				  break;
			  }
		  }
		  stream.clear();
	  }
    objects.push_back(
        new Triangle(this, v0_id - 1, v1_id - 1, v2_id - 1, 0, material_id, Arbitrary_transformation(arbitrary_transformation)));
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
    const Vector3& center_of_sphere = vertex_data[center - 1].get_vertex_position();

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
          arbitrary_transformation = scaling_transformations[index].get_transformation_matrix() * arbitrary_transformation;
          break;
        case 't':
          arbitrary_transformation = translation_transformations[index].get_transformation_matrix() * arbitrary_transformation;
          break;
        case 'r':
          arbitrary_transformation = rotation_transformations[index].get_transformation_matrix() * arbitrary_transformation;
          break;
        }
      }
      stream.clear();
    }
    objects.push_back(new Sphere(center_of_sphere, radius, material_id, Arbitrary_transformation(arbitrary_transformation)));
    element = element->NextSiblingElement("Sphere");
  }
  stream.clear();
  debug("Spheres are parsed");

  bvh = BVH::create_bvh(objects);
  //Finalize surface normals
  for(Vertex& vertex: vertex_data) {
    if(vertex.has_vertex_normal())
    {
      vertex.finalize_normal();
    }
  }
}

Scene::~Scene() {
	delete bvh;
	size_t size = meshes.size();
	for (size_t i = 0; i < size; i++) {
		delete meshes[i];
	}
}
