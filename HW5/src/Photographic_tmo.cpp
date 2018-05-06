#include "Photographic_tmo.h"

Photographic_tmo::Photographic_tmo(float image_key, float saturation_percentage,
                                   float saturation, float gamma)
    : image_key_(image_key),
      saturation_percentage_(saturation_percentage),
      saturation_(saturation),
      gamma_(gamma) {}

void Photographic_tmo::apply_tmo(const std::vector<Vector3>& input,
                                 std::vector<Vector3>& output) const {}
