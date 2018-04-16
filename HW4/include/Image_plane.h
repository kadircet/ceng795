#pragma once
#ifndef IMAGE_PLANE_H_
#define IMAGE_PLANE_H_
struct Image_plane {
  float left, right, bottom, top;
  float distance;
  int width, height;
  Image_plane(float l, float r, float b, float t, float d, int image_width,
              int image_height)
      : left(l),
        right(r),
        bottom(b),
        top(t),
        distance(d),
        width(image_width),
        height(image_height) {}
};
#endif
