#pragma once
#ifndef VECTOR3I_H_
#define VECTOR3I_H_
#include <cmath>
#include <iostream>
class Vector3i {
 public:
  float x, y, z;
  inline Vector3i() : x(0), y(0), z(0) {}
  inline Vector3i(int val) : x(val), y(val), z(val) {}
  inline Vector3i(int X, int Y, int Z) : x(X), y(Y), z(Z) {}
  //
  inline Vector3i operator+(const Vector3i& rhs) const {
    return Vector3i(x + rhs.x, y + rhs.y, z + rhs.z);
  }
  inline Vector3i operator-(const Vector3i& rhs) const {
    return Vector3i(x - rhs.x, y - rhs.y, z - rhs.z);
  }
  inline Vector3i operator*(const Vector3i& rhs) const {
    return Vector3i(x * rhs.x, y * rhs.y, z * rhs.z);
  }
  inline Vector3i operator/(const Vector3i& rhs) const {
    return Vector3i(x / rhs.x, y / rhs.y, z / rhs.z);
  }
  //
  inline Vector3i& operator+=(const Vector3i& rhs) {
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    return *this;
  }
  inline Vector3i& operator-=(const Vector3i& rhs) {
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
    return *this;
  }
  inline Vector3i& operator*=(const Vector3i& rhs) {
    x *= rhs.x;
    y *= rhs.y;
    z *= rhs.z;
    return *this;
  }
  inline Vector3i& operator/=(const Vector3i& rhs) {
    x /= rhs.x;
    y /= rhs.y;
    z /= rhs.z;
    return *this;
  }
  //
  inline Vector3i operator+(float rhs) const {
    return Vector3i(x + rhs, y + rhs, z + rhs);
  }
  inline Vector3i operator-(float rhs) const {
    return Vector3i(x - rhs, y - rhs, z - rhs);
  }
  inline Vector3i operator*(float rhs) const {
    return Vector3i(x * rhs, y * rhs, z * rhs);
  }
  inline Vector3i operator/(float rhs) const {
    return Vector3i(x / rhs, y / rhs, z / rhs);
  }
  //
  inline Vector3i& operator+=(float rhs) {
    x += rhs;
    y += rhs;
    z += rhs;
    return *this;
  }
  inline Vector3i& operator-=(float rhs) {
    x -= rhs;
    y -= rhs;
    z -= rhs;
    return *this;
  }
  inline Vector3i& operator*=(float rhs) {
    x *= rhs;
    y *= rhs;
    z *= rhs;
    return *this;
  }
  inline Vector3i& operator/=(float rhs) {
    x /= rhs;
    y /= rhs;
    z /= rhs;
    return *this;
  }
  //
  inline bool operator==(const Vector3i& rhs) const {
    return x == rhs.x && y == rhs.y && z == rhs.z;
  }
  inline bool operator!=(const Vector3i& rhs) const { return !(*this == rhs); }
  inline Vector3i operator-() const { return Vector3i(-x, -y, -z); }
  inline const float& operator[](const int i) const {
    return i == 0 ? this->x : (i == 1 ? this->y : this->z);
  }
  inline float& operator[](const int i) {
    return i == 0 ? this->x : (i == 1 ? this->y : this->z);
  }

  inline float dot(const Vector3i& rhs) const {
    return x * rhs.x + y * rhs.y + z * rhs.z;
  }
  inline Vector3i cross(const Vector3i& rhs) const {
    return Vector3i(y * rhs.z - z * rhs.y, z * rhs.x - x * rhs.z,
                    x * rhs.y - y * rhs.x);
  }
  inline float length() const { return std::sqrt(x * x + y * y + z * z); }
  inline Vector3i normalize() const { return (*this) / this->length(); }
  friend std::ostream& operator<<(std::ostream& os, const Vector3i& t);
};
inline Vector3i operator+(float a, const Vector3i& b) { return b + a; }
inline Vector3i operator-(float a, const Vector3i& b) {
  return Vector3i(a) - b;
}
inline Vector3i operator*(float a, const Vector3i& b) { return b * a; }
inline Vector3i operator/(float a, const Vector3i& b) {
  return Vector3i(a) / b;
}
#endif
