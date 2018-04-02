#pragma once
#ifndef AREA_LIGHT_H
#define AREA_LIGHT_LIGHT_H
#include "Vector3.h"
struct Area_light {
  Vector3 position;
  Vector3 intensity;
  Vector3 edge_vector_1;
  Vector3 edge_vector_2;
};
#endif
