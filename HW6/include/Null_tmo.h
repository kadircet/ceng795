#pragma once
#ifndef NULL_TMO_H_
#define NULL_TMO_H_
#include "Tonemapping_operator.h"
class Null_tmo : public Tonemapping_operator {
 public:
  void apply_tmo(const std::vector<Vector3>& input,
                 std::vector<Vector3>& output) const override;
};
#endif
