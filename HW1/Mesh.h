#pragma once
#ifndef MESH_H_
#define MESH_H_
#include <vector>
#include "Triangle.h"
struct Mesh {
  int material_id;
  int texture_id;

  std::vector<Triangle> triangles;
};
#endif
