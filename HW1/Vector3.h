#pragma once
#ifndef VECTOR3_H_
#define VECTOR3_H_
#include <cmath>
#include <iostream>
class Vector3 {
 public:
  float x, y, z;
  inline Vector3() : x(0.0f), y(0.0f), z(0.0f) {}
  inline Vector3(float val) : x(val), y(val), z(val) {}
  inline Vector3(float X, float Y, float Z) : x(X), y(Y), z(Z) {}
  //
  inline Vector3 operator+(const Vector3& rhs) const {
    return Vector3(x + rhs.x, y + rhs.y, z + rhs.z);
  }
  inline Vector3 operator-(const Vector3& rhs) const {
    return Vector3(x - rhs.x, y - rhs.y, z - rhs.z);
  }
  inline Vector3 operator*(const Vector3& rhs) const {
    return Vector3(x * rhs.x, y * rhs.y, z * rhs.z);
  }
  inline Vector3 operator/(const Vector3& rhs) const {
    return Vector3(x / rhs.x, y / rhs.y, z / rhs.z);
  }
  //
  inline Vector3& operator+=(const Vector3& rhs) {
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    return *this;
  }
  inline Vector3& operator-=(const Vector3& rhs) {
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
    return *this;
  }
  inline Vector3& operator*=(const Vector3& rhs) {
    x *= rhs.x;
    y *= rhs.y;
    z *= rhs.z;
    return *this;
  }
  inline Vector3& operator/=(const Vector3& rhs) {
    x /= rhs.x;
    y /= rhs.y;
    z /= rhs.z;
    return *this;
  }
  //
  inline Vector3 operator+(float rhs) const {
    return Vector3(x + rhs, y + rhs, z + rhs);
  }
  inline Vector3 operator-(float rhs) const {
    return Vector3(x - rhs, y - rhs, z - rhs);
  }
  inline Vector3 operator*(float rhs) const {
    return Vector3(x * rhs, y * rhs, z * rhs);
  }
  inline Vector3 operator/(float rhs) const {
    return Vector3(x / rhs, y / rhs, z / rhs);
  }
  //
  inline Vector3& operator+=(float rhs) {
    x += rhs;
    y += rhs;
    z += rhs;
    return *this;
  }
  inline Vector3& operator-=(float rhs) {
    x -= rhs;
    y -= rhs;
    z -= rhs;
    return *this;
  }
  inline Vector3& operator*=(float rhs) {
    x *= rhs;
    y *= rhs;
    z *= rhs;
    return *this;
  }
  inline Vector3& operator/=(float rhs) {
    x /= rhs;
    y /= rhs;
    z /= rhs;
    return *this;
  }
  //
  inline bool operator==(const Vector3& rhs) const {
    return x == rhs.x && y == rhs.y && z == rhs.z;
  }
  inline bool operator!=(const Vector3& rhs) const { return !(*this == rhs); }
  inline Vector3 operator-() const { return Vector3(-x, -y, -z); }
  inline const float& operator[](const int i) const {
    return i == 0 ? this->x : (i == 1 ? this->y : this->z);
  }
  inline float& operator[](const int i) {
    return i == 0 ? this->x : (i == 1 ? this->y : this->z);
  }

  inline float dot(const Vector3& rhs) const {
    return x * rhs.x + y * rhs.y + z * rhs.z;
  }
  inline Vector3 cross(const Vector3& rhs) const {
    return Vector3(y * rhs.z - z * rhs.y, z * rhs.x - x * rhs.z,
                   x * rhs.y - y * rhs.x);
  }
  inline float length() const { return std::sqrt(x * x + y * y + z * z); }
  inline Vector3 normalize() const { return (*this) / this->length(); }
  friend std::ostream& operator<<(std::ostream& os, const Vector3& t);
};
inline Vector3 operator+(float a, const Vector3& b) { return b + a; }
inline Vector3 operator-(float a, const Vector3& b) { return Vector3(a) - b; }
inline Vector3 operator*(float a, const Vector3& b) { return b * a; }
inline Vector3 operator/(float a, const Vector3& b) { return Vector3(a) / b; }
#endif
