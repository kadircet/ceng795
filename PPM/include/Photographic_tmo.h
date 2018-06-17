#pragma once
#ifndef PHOTOGRAPHIC_TMO_H_
#define PHOTOGRAPHIC_TMO_H_
#include "Tonemapping_operator.h"
class Photographic_tmo : public Tonemapping_operator {
 public:
  Photographic_tmo(float image_key = 0.18f, float saturation_percentage = 1.0f,
                   float saturation = 1.0f);
  void apply_tmo(const std::vector<Vector3>& input,
                 std::vector<Vector3>& output) const override;

 private:
  float image_key_;
  float saturation_percentage_;
  float saturation_;
  float gamma_;
};
#endif
