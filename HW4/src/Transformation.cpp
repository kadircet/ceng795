#include "Transformation.h"
#include "Vector3.h"
#include <cmath>

//Scaling
Scaling::Scaling(float x, float y, float z) {
	transformation_[0][0] = x;
	transformation_[1][1] = y;
	transformation_[2][2] = z;
	transformation_[3][3] = 1.0f;
	inverse_transformation_[0][0] = 1.0f / x;
	inverse_transformation_[1][1] = 1.0f / y;
	inverse_transformation_[2][2] = 1.0f / z;
	inverse_transformation_[3][3] = 1.0f;

	normal_transformation_ = inverse_transformation_.transpose();
}
//

//Translation 
Translation::Translation(float x, float y, float z) {
	transformation_.make_identity();
	transformation_[0][3] = x;
	transformation_[1][3] = y;
	transformation_[2][3] = z;
	inverse_transformation_.make_identity();
	inverse_transformation_[0][3] = -x;
	inverse_transformation_[1][3] = -y;
	inverse_transformation_[2][3] = -z;
	normal_transformation_ = inverse_transformation_.transpose();
}
//

//Rotation
Rotation::Rotation(float angle, float x, float y, float z) {
	const Vector3 u = Vector3(x, y, z).normalize();
	const Vector3 v = ((x != 0.0f || y != 0.0f) ? Vector3(-u.y, u.x, 0.0f) : Vector3(0.0f, 1.0f, 0.0f))
		.normalize();
	const Vector3 w = u.cross(v);

	Matrix4x4 m;
	m[0][0] = u.x;
	m[0][1] = u.y;
	m[0][2] = u.z;
	m[1][0] = v.x;
	m[1][1] = v.y;
	m[1][2] = v.z;
	m[2][0] = w.x;
	m[2][1] = w.y;
	m[2][2] = w.z;
	m[3][3] = 1.0f;

	Matrix4x4 rot;
	rot[0][0] = 1.0f;
	rot[1][1] = std::cos(angle);
	rot[1][2] = -std::sin(angle);
	rot[2][1] = -rot[1][2];
	rot[2][2] = rot[1][1];
	rot[3][3] = 1.0f;
	Matrix4x4 inv_rot;
	inv_rot[0][0] = 1.0f;
	inv_rot[1][1] = std::cos(-angle);
	inv_rot[1][2] = -std::sin(-angle);
	inv_rot[2][1] = -inv_rot[1][2];
	inv_rot[2][2] = inv_rot[1][1];
	inv_rot[3][3] = 1.0f;
	Matrix4x4 m_transpose = m.transpose();
	transformation_ = m_transpose * (rot * m);
	inverse_transformation_ = m_transpose * (inv_rot * m);

	normal_transformation_ = inverse_transformation_.transpose();
}
//

//ArbitraryTransformation
Arbitrary_transformation::Arbitrary_transformation(const Matrix4x4& transformation_matrix) {
	transformation_ = transformation_matrix;
	if (!(transformation_matrix.invert_matrix(inverse_transformation_)))
	{
		std::cerr << "NOT INVERTIBLE MATRIX!" << std::endl;
		exit(-1);
	}
	normal_transformation_= inverse_transformation_.transpose();
}
