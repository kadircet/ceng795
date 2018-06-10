#include "Photographic_tmo.h"
#include <algorithm>

Photographic_tmo::Photographic_tmo(float image_key, float saturation_percentage,
                                   float saturation)
    : image_key_(image_key),
      saturation_percentage_(saturation_percentage),
      saturation_(saturation) {}

void Photographic_tmo::apply_tmo(const std::vector<Vector3>& input,
                                 std::vector<Vector3>& output) const {
  float log_average_luminance = 0;
  size_t size = input.size();
  std::vector<float> luminances;
  const float epsilon = 0.001;
  for (size_t i = 0; i < size; i++) {
    const Vector3& color = input[i];
    float lum = (0.21f * color.x + 0.72f * color.y + 0.07f * color.z) / 255.0f;
    luminances.push_back(lum);
    log_average_luminance += std::log(lum + epsilon);
  }
  log_average_luminance = std::exp(log_average_luminance / (float)size);
  for (size_t i = 0; i < size; i++) {
    luminances[i] = image_key_ * luminances[i] / log_average_luminance;
  }
  std::vector<float> sorted_luminances = luminances;
  std::sort(sorted_luminances.begin(), sorted_luminances.end());
  float white_count = (float)size * saturation_percentage_ / 100.0f;
  float white_lum = sorted_luminances[size - (size_t)white_count];
  for (size_t i = 0; i < size; i++) {
    luminances[i] =
        (luminances[i] * (1 + luminances[i] / (white_lum * white_lum))) /
        (luminances[i] + 1);
    const Vector3& color = input[i];
    float lum_i = (0.21f * color.x + 0.72f * color.y + 0.07f * color.z);
    float r_f = std::pow((color.x / lum_i), saturation_) * luminances[i];
    float g_f = std::pow((color.y / lum_i), saturation_) * luminances[i];
    float b_f = std::pow((color.z / lum_i), saturation_) * luminances[i];
    float r_d = std::min(1.0f, std::max(0.0f, r_f)) * 255.0f;
    float g_d = std::min(1.0f, std::max(0.0f, g_f)) * 255.0f;
    float b_d = std::min(1.0f, std::max(0.0f, b_f)) * 255.0f;
    output.push_back(Vector3(r_d, g_d, b_d));
  }
}
