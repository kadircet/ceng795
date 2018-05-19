#include "Null_tmo.h"
#include <algorithm>
#include "Vector3.h"

void Null_tmo::apply_tmo(const std::vector<Vector3>& input,
                         std::vector<Vector3>& output) const {
  for (const Vector3& color : input) {
    output.push_back(Vector3(std::min(255.0f, std::max(0.0f, color.x)),
                             std::min(255.0f, std::max(0.0f, color.y)),
                             std::min(255.0f, std::max(0.0f, color.z))));
  }
}
