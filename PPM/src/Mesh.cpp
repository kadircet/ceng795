#include "Mesh.h"
#include <sstream>
#include "Mesh_triangle.h"
#include "Scene.h"
void Mesh::load_meshes_from_xml(
    const Scene* scene, tinyxml2::XMLElement* element,
    std::vector<Mesh*>& meshes, std::vector<Vertex>& vertex_data,
    std::vector<Scaling>& scaling_transformations,
    std::vector<Translation>& translation_transformations,
    std::vector<Rotation>& rotation_transformations) {
  element = element->FirstChildElement("Mesh");
  std::stringstream stream;
  while (element) {
    auto child = element->FirstChildElement("Material");
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
    const char* binary_file = child->Attribute("binaryFile");
    if (ply_file) {
      // parse_ply_tinyply(std::string(ply_file), vertex_data, triangles,
      //                  vertex_offset, texture_offset, material_id,
      //                  texture_id, triangle_shading_mode);
    } else if (binary_file) {
      // vertex_offset = child->IntAttribute("vertexOffset", 0);
      // parse_binary_facedata(std::string(binary_file), triangles,
      // vertex_offset,
      //                      texture_offset, material_id, texture_id,
      //                      triangle_shading_mode, zero_based_indexing);
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
            scene, v0_id, v1_id, v2_id, vertex_offset, texture_offset,
            material_id, texture_id, triangle_shading_mode);
        float area = triangle->get_surface_area();
        const Vector3& surface_normal = triangle->normal;
        vertex_data[triangle->vertex_index_0 + triangle->vertex_offset]
            .add_vertex_normal(surface_normal, area);
        vertex_data[triangle->vertex_index_1 + triangle->vertex_offset]
            .add_vertex_normal(surface_normal, area);
        vertex_data[triangle->vertex_index_2 + triangle->vertex_offset]
            .add_vertex_normal(surface_normal, area);
        triangles.push_back(triangle);
      }
    }
    stream.clear();
    meshes.push_back(
        new Mesh(material_id, texture_id, triangles,
                 Arbitrary_transformation(arbitrary_transformation)));
    element = element->NextSiblingElement("Mesh");
  }
  stream.clear();
  std::cout << "Meshes are parsed" << std::endl;
}

void Mesh_instance::create_mesh_instances_for_meshes(
    std::vector<Mesh*>& meshes, std::vector<Shape*>& objects,
    std::vector<Material>& materials) {
  for (Mesh* mesh : meshes) {
    const Material& material = materials[mesh->get_material_id()];
    bool is_refractive = material.material_type == mt_refractive;
    objects.push_back(new Mesh_instance(mesh->get_material_id(),
                                        mesh->texture_id, mesh,
                                        mesh->base_transform, is_refractive));
  }
  std::cout << "Base mesh instances are created" << std::endl;
}
void Mesh_instance::load_mesh_instances_from_xml(
    tinyxml2::XMLElement* element, std::vector<Mesh*>& meshes,
    std::vector<Shape*>& objects, std::vector<Material>& materials,
    std::vector<Scaling>& scaling_transformations,
    std::vector<Translation>& translation_transformations,
    std::vector<Rotation>& rotation_transformations) {
  element = element->FirstChildElement("MeshInstance");
  std::stringstream stream;

  while (element) {
    Mesh* base_mesh = meshes[element->IntAttribute("baseMeshId") - 1];
    auto child = element->FirstChildElement("Material");
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

    // TODO: Until finding an elegant way, different textures for different
    // mesh instances are not supported. Since bump_map and perlin_noise
    // calculations are done in mesh_triangle of base_mesh Normal textures may
    // work
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
    const Material& material = materials[material_id];
    bool is_refractive = material.material_type == mt_refractive;
    objects.push_back(new Mesh_instance(
        material_id, texture_id, base_mesh,
        Arbitrary_transformation(arbitrary_transformation), is_refractive));
    element = element->NextSiblingElement("MeshInstance");
  }
  stream.clear();
  std::cout << "MeshInstances are parsed" << std::endl;
}
