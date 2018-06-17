#include "Vector3.h"

Vector3 Vector3::zero_vector = Vector3(0.0f);
std::ostream& operator<<(std::ostream& os, const Vector3& v) {
  os << "Vector3(" << v.x << ", " << v.y << ", " << v.z << ")";
  return os;
}
