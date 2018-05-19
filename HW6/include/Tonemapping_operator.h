#pragma once
#ifndef TONEMAPPING_OPERATOR_H_
#define TONEMAPPING_OPERATOR_H_
#include <vector>
#include "Vector3.h"
class Tonemapping_operator {
 public:
  virtual void apply_tmo(const std::vector<Vector3>& input,
                         std::vector<Vector3>& output) const = 0;
};
#endif
