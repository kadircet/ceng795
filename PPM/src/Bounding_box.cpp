#include "Bounding_box.h"
#include <algorithm>
#include "Matrix4x4.h"
#include "Transformation.h"

void Bounding_box::expand(const Bounding_box& bounding_box) {
  min_corner.x = std::min(min_corner.x, bounding_box.min_corner.x);
  min_corner.y = std::min(min_corner.y, bounding_box.min_corner.y);
  min_corner.z = std::min(min_corner.z, bounding_box.min_corner.z);

  max_corner.x = std::max(max_corner.x, bounding_box.max_corner.x);
  max_corner.y = std::max(max_corner.y, bounding_box.max_corner.y);
  max_corner.z = std::max(max_corner.z, bounding_box.max_corner.z);

  delta = max_corner - min_corner;
  center = (max_corner + min_corner) / 2.;
}

void Bounding_box::fit(const Vector3& point) {
  min_corner.x = std::min(min_corner.x, point.x);
  min_corner.y = std::min(min_corner.y, point.y);
  min_corner.z = std::min(min_corner.z, point.z);

  max_corner.x = std::max(max_corner.x, point.x);
  max_corner.y = std::max(max_corner.y, point.y);
  max_corner.z = std::max(max_corner.z, point.z);

  delta = max_corner - min_corner;
  center = (max_corner + min_corner) / 2.;
}

float Bounding_box::intersect(const Ray& ray) const {
  float tmin = -kInf;
  float tmax = kInf;
  const Vector3& origin = ray.o;
  const Vector3& direction = ray.d;
  for (int i = 0; i < 3; i++) {
    if (fabs(direction[i]) < intersection_test_epsilon) continue;
    float ti_min = (min_corner[i] - origin[i]) / direction[i];
    float ti_max = (max_corner[i] - origin[i]) / direction[i];
    if (direction[i] < 0) std::swap(ti_min, ti_max);
    if (ti_min > tmin) tmin = ti_min;
    if (ti_max < tmax) tmax = ti_max;
    if (tmin > tmax) return kInf;
  }

  if (tmin > 0.0f) {
    return tmin;
  } else {
    return tmax;
  }
}

Bounding_box Bounding_box::apply_transform(const Bounding_box& bounding_box,
                                           const Transformation& transform) {
  const Matrix4x4& transformation_matrix =
      transform.get_transformation_matrix();
  const Vector3& min = bounding_box.min_corner;
  const Vector3& max = bounding_box.max_corner;

  const Vector3 v000 =
      transformation_matrix.multiply(Vector3(min.x, min.y, min.z));
  const Vector3 v001 =
      transformation_matrix.multiply(Vector3(min.x, min.y, max.z));
  const Vector3 v010 =
      transformation_matrix.multiply(Vector3(min.x, max.y, min.z));
  const Vector3 v011 =
      transformation_matrix.multiply(Vector3(min.x, max.y, max.z));
  const Vector3 v100 =
      transformation_matrix.multiply(Vector3(max.x, min.y, min.z));
  const Vector3 v101 =
      transformation_matrix.multiply(Vector3(max.x, min.y, max.z));
  const Vector3 v110 =
      transformation_matrix.multiply(Vector3(max.x, max.y, min.z));
  const Vector3 v111 =
      transformation_matrix.multiply(Vector3(max.x, max.y, max.z));

  /// v000
  Vector3 min_c = v000;
  Vector3 max_c = v000;
  /// v001
  min_c.x = std::min(min_c.x, v001.x);
  min_c.y = std::min(min_c.y, v001.y);
  min_c.z = std::min(min_c.z, v001.z);

  max_c.x = std::max(max_c.x, v001.x);
  max_c.y = std::max(max_c.y, v001.y);
  max_c.z = std::max(max_c.z, v001.z);
  /// v010
  min_c.x = std::min(min_c.x, v010.x);
  min_c.y = std::min(min_c.y, v010.y);
  min_c.z = std::min(min_c.z, v010.z);

  max_c.x = std::max(max_c.x, v010.x);
  max_c.y = std::max(max_c.y, v010.y);
  max_c.z = std::max(max_c.z, v010.z);
  /// v011
  min_c.x = std::min(min_c.x, v011.x);
  min_c.y = std::min(min_c.y, v011.y);
  min_c.z = std::min(min_c.z, v011.z);

  max_c.x = std::max(max_c.x, v011.x);
  max_c.y = std::max(max_c.y, v011.y);
  max_c.z = std::max(max_c.z, v011.z);
  /// v100
  min_c.x = std::min(min_c.x, v100.x);
  min_c.y = std::min(min_c.y, v100.y);
  min_c.z = std::min(min_c.z, v100.z);

  max_c.x = std::max(max_c.x, v100.x);
  max_c.y = std::max(max_c.y, v100.y);
  max_c.z = std::max(max_c.z, v100.z);
  /// v101
  min_c.x = std::min(min_c.x, v101.x);
  min_c.y = std::min(min_c.y, v101.y);
  min_c.z = std::min(min_c.z, v101.z);

  max_c.x = std::max(max_c.x, v101.x);
  max_c.y = std::max(max_c.y, v101.y);
  max_c.z = std::max(max_c.z, v101.z);
  /// v110
  min_c.x = std::min(min_c.x, v110.x);
  min_c.y = std::min(min_c.y, v110.y);
  min_c.z = std::min(min_c.z, v110.z);

  max_c.x = std::max(max_c.x, v110.x);
  max_c.y = std::max(max_c.y, v110.y);
  max_c.z = std::max(max_c.z, v110.z);
  /// v111
  min_c.x = std::min(min_c.x, v111.x);
  min_c.y = std::min(min_c.y, v111.y);
  min_c.z = std::min(min_c.z, v111.z);

  max_c.x = std::max(max_c.x, v111.x);
  max_c.y = std::max(max_c.y, v111.y);
  max_c.z = std::max(max_c.z, v111.z);
  return Bounding_box(min_c, max_c);
}
