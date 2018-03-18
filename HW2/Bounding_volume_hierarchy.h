#ifndef BOUNDING_VOLUME_HIERARCHY_H_
#define BOUNDING_VOLUME_HIERARCHY_H_
#include <vector>
#include "Bounding_box.h"
#include "Shape.h"
#include "Vector3.h"
class BVH : public Shape {
 public:
  static Shape* create_bvh(std::vector<Shape*>& objects) {
    int size = objects.size();
    if (size == 0) {
      return NULL;
    } else if (size == 1) {
      return objects[0];
    } else {
      return new BVH(objects, 0, size, 0);
    }
  }
  BVH(std::vector<Shape*>& objects, int start, int end, int dimension);
  ~BVH() {
    if (left) {
      delete left;
    }
    if (right) {
      delete right;
    }
  }
  Hit_data intersect(const Ray& ray) const override;
  int get_material_id() const override { return -1; }
  const Bounding_box& get_bounding_box() const override { return bounding_box; }
  void print_debug(int indentation) const override {
    for (int index = 0; index < indentation; index++) {
      std::cout << "\t";
    }
    std::cout << "BVH:" << std::endl;
    left->print_debug(indentation + 1);
    right->print_debug(indentation + 1);
  }
  Shape* left;
  Shape* right;
  Bounding_box bounding_box;
};
#endif
