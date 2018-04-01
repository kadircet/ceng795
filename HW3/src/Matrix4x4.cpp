#include "Matrix4x4.h"
#include "Vector3.h"

Matrix4x4::Matrix4x4(bool is_identity) {
	if (is_identity) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				elements_[i][j] = ((i == j) ? 1.0f : 0.0f);
			}
		}
	}
	else {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				elements_[i][j] = 0.0f;
			}
		}	
	}
}

Matrix4x4::Matrix4x4(float element) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			elements_[i][j] = element;
		}
	}
}

Vector3 Matrix4x4::multiply(const Vector3& rhs, bool is_vector) const {
	Vector3 result;
	if (is_vector) {
		for (int i = 0; i < 3; i++) {
			result[i] = 0;
			for (int j = 0; j < 3; j++) {
				result[i] += elements_[i][j] * rhs[j];
			}
		}
	}
	else {
		//point
		for (int i = 0; i < 3; i++) {
			result[i] = elements_[i][3];
			for (int j = 0; j < 3; j++) {
				result[i] += elements_[i][j] * rhs[j];
			}
		}
	}
	return result;
}

Matrix4x4 Matrix4x4::operator*(const Matrix4x4& rhs) const {
	Matrix4x4 result;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				result[i][j] += elements_[i][k] * rhs[k][j];
			}
		}
	}
	return result;
}

Matrix4x4 Matrix4x4::operator*=(const Matrix4x4& rhs) {
	float result[4][4];
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			result[i][j] = 0;
			for (int k = 0; k < 4; k++) {
				result[i][j] += elements_[i][k] * rhs[k][j];
			}
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			elements_[i][j] = result[i][j];
		}
	}
	return *this;
}

Matrix4x4 Matrix4x4::transpose() const {
	Matrix4x4 t;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			t[i][j] = elements_[j][i];
		}
	}
	return t;
}

void Matrix4x4::make_identity() {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			elements_[i][j] = ((i == j) ? 1.0f : 0.0f);
		}
	}
}

bool Matrix4x4::invert_matrix(Matrix4x4& inverse_out) const
{
	float inv[16], m[16];
	float det;
	int i;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			m[i * 4 + j] = elements_[i][j];
		}
	}
	inv[0] = m[5] * m[10] * m[15] -
		m[5] * m[11] * m[14] -
		m[9] * m[6] * m[15] +
		m[9] * m[7] * m[14] +
		m[13] * m[6] * m[11] -
		m[13] * m[7] * m[10];

	inv[4] = -m[4] * m[10] * m[15] +
		m[4] * m[11] * m[14] +
		m[8] * m[6] * m[15] -
		m[8] * m[7] * m[14] -
		m[12] * m[6] * m[11] +
		m[12] * m[7] * m[10];

	inv[8] = m[4] * m[9] * m[15] -
		m[4] * m[11] * m[13] -
		m[8] * m[5] * m[15] +
		m[8] * m[7] * m[13] +
		m[12] * m[5] * m[11] -
		m[12] * m[7] * m[9];

	inv[12] = -m[4] * m[9] * m[14] +
		m[4] * m[10] * m[13] +
		m[8] * m[5] * m[14] -
		m[8] * m[6] * m[13] -
		m[12] * m[5] * m[10] +
		m[12] * m[6] * m[9];

	inv[1] = -m[1] * m[10] * m[15] +
		m[1] * m[11] * m[14] +
		m[9] * m[2] * m[15] -
		m[9] * m[3] * m[14] -
		m[13] * m[2] * m[11] +
		m[13] * m[3] * m[10];

	inv[5] = m[0] * m[10] * m[15] -
		m[0] * m[11] * m[14] -
		m[8] * m[2] * m[15] +
		m[8] * m[3] * m[14] +
		m[12] * m[2] * m[11] -
		m[12] * m[3] * m[10];

	inv[9] = -m[0] * m[9] * m[15] +
		m[0] * m[11] * m[13] +
		m[8] * m[1] * m[15] -
		m[8] * m[3] * m[13] -
		m[12] * m[1] * m[11] +
		m[12] * m[3] * m[9];

	inv[13] = m[0] * m[9] * m[14] -
		m[0] * m[10] * m[13] -
		m[8] * m[1] * m[14] +
		m[8] * m[2] * m[13] +
		m[12] * m[1] * m[10] -
		m[12] * m[2] * m[9];

	inv[2] = m[1] * m[6] * m[15] -
		m[1] * m[7] * m[14] -
		m[5] * m[2] * m[15] +
		m[5] * m[3] * m[14] +
		m[13] * m[2] * m[7] -
		m[13] * m[3] * m[6];

	inv[6] = -m[0] * m[6] * m[15] +
		m[0] * m[7] * m[14] +
		m[4] * m[2] * m[15] -
		m[4] * m[3] * m[14] -
		m[12] * m[2] * m[7] +
		m[12] * m[3] * m[6];

	inv[10] = m[0] * m[5] * m[15] -
		m[0] * m[7] * m[13] -
		m[4] * m[1] * m[15] +
		m[4] * m[3] * m[13] +
		m[12] * m[1] * m[7] -
		m[12] * m[3] * m[5];

	inv[14] = -m[0] * m[5] * m[14] +
		m[0] * m[6] * m[13] +
		m[4] * m[1] * m[14] -
		m[4] * m[2] * m[13] -
		m[12] * m[1] * m[6] +
		m[12] * m[2] * m[5];

	inv[3] = -m[1] * m[6] * m[11] +
		m[1] * m[7] * m[10] +
		m[5] * m[2] * m[11] -
		m[5] * m[3] * m[10] -
		m[9] * m[2] * m[7] +
		m[9] * m[3] * m[6];

	inv[7] = m[0] * m[6] * m[11] -
		m[0] * m[7] * m[10] -
		m[4] * m[2] * m[11] +
		m[4] * m[3] * m[10] +
		m[8] * m[2] * m[7] -
		m[8] * m[3] * m[6];

	inv[11] = -m[0] * m[5] * m[11] +
		m[0] * m[7] * m[9] +
		m[4] * m[1] * m[11] -
		m[4] * m[3] * m[9] -
		m[8] * m[1] * m[7] +
		m[8] * m[3] * m[5];

	inv[15] = m[0] * m[5] * m[10] -
		m[0] * m[6] * m[9] -
		m[4] * m[1] * m[10] +
		m[4] * m[2] * m[9] +
		m[8] * m[1] * m[6] -
		m[8] * m[2] * m[5];

	det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

	if (det == 0.0f)
		return false;

	det = 1.0f / det;

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			inverse_out[i][j] = inv[i * 4 + j]*det;
		}
	}

	return true;
}