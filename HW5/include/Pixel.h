#ifndef PIXEL_H_
#define PIXEL_H_
#include <mutex>
#include <thread>
#include "Vector3.h"
#include <algorithm>

class Pixel {
 public:
  Pixel() : color(0.0f), weight(0.0f) {}
  Vector3 color;
  float weight;
  void add_color(const Vector3& color, float weight) {
    std::lock_guard<std::mutex> lock(mutex_);
    this->color += (color * weight);
    this->weight += weight;
  }
  Vector3 get_color() {
    std::lock_guard<std::mutex> lock(mutex_);
    Vector3 result_color;
    if (weight == 0) {
      return Vector3(0.0f);
    } else {
      result_color = color / weight;
      return Vector3(std::min(255.0f, std::max(0.0f, result_color.x)),
					std::min(255.0f, std::max(0.0f, result_color.y)),
					std::min(255.0f, std::max(0.0f, result_color.z)));
    }
  }

 private:
  std::mutex mutex_;
};
#endif
