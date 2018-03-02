#include "Vector3i.h"

std::ostream& operator<<(std::ostream& os, const Vector3i& v) {
  os << "Vector3i(" << v.x << ", " << v.y << ", " << v.z << ")";
  return os;
}
