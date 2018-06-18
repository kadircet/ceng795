#pragma once
#ifndef HIT_POINT_H_
#define HIT_POINT_H_
#include <mutex>
#include "Material.h"
#include "Vector3.h"
class Hit_point {
 public:
  // Color calculation
  Material material;
  Vector3 attenuation;
  Vector3 w_o;
  Vector3 normal;
  //
  Vector3 position;
  Vector3 flux;
  float radius_squared;
  unsigned int n;
  int pixel;
  float pixel_weight;

  std::mutex mutex;
};
#endif
