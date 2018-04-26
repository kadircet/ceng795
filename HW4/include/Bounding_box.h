#ifndef BOUNDING_BOX_H_
#define BOUNDING_BOX_H_
#include "Ray.h"
#include "Vector3.h"
class Transformation;
class Bounding_box {
 public:
  static Bounding_box apply_transform(const Bounding_box& bounding_box,
                                      const Transformation& transform);
  Bounding_box()
      : min_corner(kInf, kInf, kInf), max_corner(-kInf, -kInf, -kInf) {}
  Bounding_box(const Vector3& min, const Vector3& max)
      : min_corner(min),
        max_corner(max),
        delta(max - min),
        center((max + min) / 2) {}
  void expand(const Bounding_box& bounding_box);
  float intersect(const Ray& ray) const;
  Vector3 min_corner;
  Vector3 max_corner;
  Vector3 delta;
  Vector3 center;
  int max_dimension() const {
    if (delta.x > delta.y) {
      if (delta.x > delta.z) return 0;
      return 2;
    }
    if (delta.y > delta.z) return 1;
    return 2;
  }
};
#endif
