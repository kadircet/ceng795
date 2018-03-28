#include <chrono>
#include <iostream>
#include <thread>
#include "Pixel.h"
#include "Scene.h"
#include "lodepng.h"
#include "timeutil.h"
#define THREAD_MULTIPLIER 1

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Please provide scene file as argument" << std::endl;
    return 1;
  }
  Scene scene(argv[1]);
  std::cout << "Scene is parsed" << std::endl;
  const int thread_count =
      std::thread::hardware_concurrency() * THREAD_MULTIPLIER;
  const int camera_count = (int) scene.cameras.size();
  for (int index = 0; index < camera_count; index++) {
    const Camera& camera = scene.cameras[index];
    const Image_plane& image_plane = camera.get_image_plane();
    const int width = image_plane.width;
    const int height = image_plane.height;
    Pixel* pixels = new Pixel[width * height];
    auto start = std::chrono::system_clock::now();
    if (thread_count == 0 || height < thread_count) {
      std::cout << "Starting rendering on #1 thread(s)" << std::endl;
      scene.render_image(index, pixels, 0);
    } else {
      std::cout << "Starting rendering on #" << thread_count << " thread(s)"
                << std::endl;
      std::thread* threads = new std::thread[thread_count];
      for (int i = 0; i < thread_count; i++) {
        threads[i] = std::thread(&Scene::render_image, &scene, index, pixels, i,
                                 thread_count);
      }
      for (int i = 0; i < thread_count; i++) threads[i].join();
      delete[] threads;
    }
    auto end = std::chrono::system_clock::now();
    // Encode the image
    unsigned char* image = new unsigned char[width * height * 4];

    int idx = 0;
    for (int j = 0; j < height; j++) {
      for (int i = 0; i < width; i++) {
        const Vector3 pixel = pixels[j * width + i].get_color();
        image[idx++] = (int) pixel.x;
        image[idx++] = (int) pixel.y;
        image[idx++] = (int) pixel.z;
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
    std::cout << camera.get_filename() << "(" << width << "x" << height
              << ") is saved in: ";
    print_time_diff(std::cout, start, end);
    std::cout << std::endl;
    delete[] image;
  }
  return 0;
}
