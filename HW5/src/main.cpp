#include <chrono>
#include <iostream>
#include <opencv2/highgui/highgui.hpp>
#include <thread>
#include "Pixel.h"
#include "Scene.h"
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
  const int camera_count = (int)scene.cameras.size();
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
    std::vector<Vector3> pixel_colors;
    std::vector<Vector3> tonemapped_colors;
    const Tonemapping_operator* tmo = camera.get_tmo();
    for (int j = 0; j < height; j++) {
      for (int i = 0; i < width; i++) {
        pixel_colors.push_back(pixels[j * width + i].get_color());
      }
    }
    tmo->apply_tmo(pixel_colors, tonemapped_colors);
    // Encode the image
    cv::Mat image = cv::Mat(height, width, CV_32FC3);

    int idx = 0;
    for (int j = 0; j < height; j++) {
      for (int i = 0; i < width; i++) {
        const Vector3 pixel = tonemapped_colors[j * width + i];
        image.at<cv::Vec3f>(j, i)[0] = pixel.z;
        image.at<cv::Vec3f>(j, i)[1] = pixel.y;
        image.at<cv::Vec3f>(j, i)[2] = pixel.x;
      }
    }
    delete[] pixels;
    cv::imwrite(camera.get_filename(), image);
    std::cout << camera.get_filename() << "(" << width << "x" << height
              << ") is saved in: ";
    print_time_diff(std::cout, start, end);
    std::cout << std::endl;
  }
  return 0;
}
