#include "Material.h"
#include <sstream>
void Material::load_materials_from_xml(tinyxml2::XMLElement* element,
                                       std::vector<Material>& materials) {
  element = element->FirstChildElement("Material");
  std::stringstream stream;
  Material material;
  while (element) {
    auto child = element->FirstChildElement("DiffuseReflectance");
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
    material.brdf_id = element->IntAttribute("BRDF", 0) - 1;

    stream >> material.diffuse.x >> material.diffuse.y >> material.diffuse.z;
    stream >> material.specular.x >> material.specular.y >> material.specular.z;
    stream >> material.mirror.x >> material.mirror.y >> material.mirror.z;
    stream >> material.phong_exponent;
    stream >> material.transparency.x >> material.transparency.y >>
        material.transparency.z;
    stream >> material.refraction_index;

    if (material.mirror != Vector3::zero_vector) {
      material.material_type = mt_mirror;
    } else if (material.transparency != Vector3::zero_vector) {
      material.material_type = mt_refractive;
    } else {
      material.material_type = mt_diffuse;
    }

    if (element->BoolAttribute("degamma", false)) {
      material.diffuse.x = pow(material.diffuse.x, 2.2f);
      material.diffuse.y = pow(material.diffuse.y, 2.2f);
      material.diffuse.z = pow(material.diffuse.z, 2.2f);
      material.specular.x = pow(material.specular.x, 2.2f);
      material.specular.y = pow(material.specular.y, 2.2f);
      material.specular.z = pow(material.specular.z, 2.2f);
    }
    materials.push_back(material);
    element = element->NextSiblingElement("Material");
  }
  stream.clear();
  std::cout << "Materials are parsed" << std::endl;
}
