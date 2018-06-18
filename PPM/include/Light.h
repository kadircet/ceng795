#ifndef LIGHT_H_
#define LIGHT_H_
#include "Ray.h"
#include "Vector3.h"
class Light {
 public:
  /*virtual Vector3 direction_and_distance(const Vector3& from_point,
                                         const Vector3& normal, float& distance,
                                         float& probability) const = 0;

  // Incoming radiance to the point from the light
  virtual Vector3 incoming_radiance(const Vector3& from_point_to_light,
                                    float probability) const = 0;*/
  virtual void generate_photon(Ray& photon_ray, Vector3& flux,
                               int photon_id) const = 0;

  static int primes[61];
  static inline int rev(const int i, const int p) {
    if (i == 0) {
      return i;
    } else {
      return p - i;
    }
  }
  static float hal(const int b, int j) {
    const int p = primes[b];
    float h = 0.0f;
    float f = 1.0f / (float)p;
    float fct = f;
    while (j > 0) {
      h += rev(j % p, p) * fct;
      j /= p;
      fct *= f;
    }
    return h;
  }
};
#endif
