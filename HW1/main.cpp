#include <iostream>
#include <thread>
#include "Scene.h"
#include "lodepng/lodepng.h"

#define THREAD_MULTIPLIER 0
int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Please provise scene file as argument" << std::endl;
    return 1;
  }
  Scene scene(argv[1]);
  std::cout << "Scene is constructed" << std::endl;
  const int thread_count =
      std::thread::hardware_concurrency() * THREAD_MULTIPLIER;
  const int camera_count = scene.cameras.size();
  for (int index = 0; index < camera_count; index++) {
    const Camera& camera = scene.cameras[index];
    const Image_plane& image_plane = camera.get_image_plane();
    const int width = image_plane.width;
    const int height = image_plane.height;
    Vector3i* pixels = new Vector3i[width * height];
    if (thread_count == 0 || height < thread_count) {
      scene.render_image(index, pixels, 0);
    } else {
      std::thread* threads = new std::thread[thread_count];
      for (int i = 0; i < thread_count; i++) {
        threads[i] = std::thread(&Scene::render_image, &scene, index, pixels, i,
                                 thread_count);
      }
      for (int i = 0; i < thread_count; i++) threads[i].join();
      delete[] threads;
    }
    // Encode the image
    unsigned char* image = new unsigned char[width * height * 4];

    int idx = 0;
    for (int j = 0; j < height; j++) {
      for (int i = 0; i < width; i++) {
        const Vector3i& pixel = pixels[j * width + i];
        image[idx++] = pixel.x;
        image[idx++] = pixel.y;
        image[idx++] = pixel.z;
        image[idx++] = 255;
      }
    }
    delete[] pixels;
    unsigned error =
        lodepng::encode(camera.get_filename().c_str(), image, width, height);

    // if there's an error, display it
    if (error)
      std::cout << "encoder error " << error << ": "
                << lodepng_error_text(error) << std::endl;
    delete[] image;
  }
  return 0;
}
