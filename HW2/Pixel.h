#ifndef PIXEL_H_
#define PIXEL_H_
#include <mutex>
#include <thread>
#include "Vector3.h"
#include "Vector3i.h"
class Pixel {
 public:
  Vector3 color;
  float weight;
  void add_color(const Vector3& color, float weight) {
    std::lock_guard<std::mutex> lock(mutex_);
    this->color += (color * weight);
    this->weight += weight;
  }
  Vector3i get_color() {
    std::lock_guard<std::mutex> lock(mutex_);
    Vector3 result_color;
    if (weight == 0) {
      return Vector3i(0);
    } else {
      result_color = color / weight;
      return Vector3i(min(255, max(0, int(result_color.x))),
                      min(255, max(0, int(result_color.y))),
                      min(255, max(0, int(result_color.z))));
    }
  }

 private:
  std::mutex mutex_;
  inline int max(int a, int b) { return a > b ? a : b; }

  inline int min(int a, int b) { return a < b ? a : b; }
};
#endif
