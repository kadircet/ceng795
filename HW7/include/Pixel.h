#ifndef PIXEL_H_
#define PIXEL_H_
#include <algorithm>
#include <mutex>
#include <thread>
#include "Vector3.h"

class Pixel {
 public:
  Pixel() : color(0.0f), weight(0.0f) {}
  Vector3 color;
  float weight;
  void add_color(const Vector3& color, float weight) {
    std::lock_guard<std::mutex> lock(mutex_);
    this->color += (color * 255.0f * weight);
    this->weight += weight;
  }
  Vector3 get_color() {
    std::lock_guard<std::mutex> lock(mutex_);
    Vector3 result_color;
    if (weight == 0) {
      return Vector3(0.0f);
    } else {
      return color / weight / 255.0f;
    }
  }

 private:
  std::mutex mutex_;
};
#endif
