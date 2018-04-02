#pragma once

#ifndef MATRIX4X4_H_
#define MATRIX4X4_H_
#include <iostream>
class Vector3;

class Matrix4x4 {
public:
	float* operator[](int index) { return elements_[index]; }
	const float* operator[](int index) const { return elements_[index]; }
	Matrix4x4(bool is_identity = false);
	Matrix4x4(float element);
	Vector3 multiply(const Vector3& rhs, bool is_vector = false) const;
	Matrix4x4 operator*(const Matrix4x4& rhs) const;
	Matrix4x4 operator*=(const Matrix4x4& rhs);
	Matrix4x4 transpose() const;
	void make_identity();
	bool invert_matrix(Matrix4x4& inverse_out) const;
  friend std::ostream& operator<<(std::ostream& os, const Matrix4x4& m);
private:
	float elements_[4][4];
};
#endif
