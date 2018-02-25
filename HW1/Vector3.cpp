#include "Vector3.h"

std::ostream &operator<<(std::ostream &os, const Vector3 &v) {
  os << "Vector3(" << v.x << ", " << v.y << ", " << v.z << ")";
  return os;
}
