#include "Bounding_box.h"
void Bounding_box::expand(const Bounding_box& bounding_box) {
  min_corner.x = fmin(min_corner.x, bounding_box.min_corner.x);
  min_corner.y = fmin(min_corner.y, bounding_box.min_corner.y);
  min_corner.z = fmin(min_corner.z, bounding_box.min_corner.z);

  max_corner.x = fmax(max_corner.x, bounding_box.max_corner.x);
  max_corner.y = fmax(max_corner.y, bounding_box.max_corner.y);
  max_corner.z = fmax(max_corner.z, bounding_box.max_corner.z);

  delta = max_corner - min_corner;
  center = (max_corner + min_corner) / 2.;
}

float Bounding_box::intersect(const Ray& ray) const {
  float tmin = -kInf;
  float tmax = kInf;
  const Vector3& origin = ray.o;
  const Vector3& direction = ray.d;
  for (int i = 0; i < 3; i++) {
    if (fabs(direction[i]) < kEpsilon) continue;
    float ti_min = (min_corner[i] - origin[i]) / direction[i];
    float ti_max = (max_corner[i] - origin[i]) / direction[i];
    if (direction[i] < 0) std::swap(ti_min, ti_max);
    if (ti_min > tmin) tmin = ti_min;
    if (ti_max < tmax) tmax = ti_max;
    if (tmin > tmax) return kInf;
  }

  if (tmin > kEpsilon) {
    return tmin;
  } else {
    return tmax;
  }
}
