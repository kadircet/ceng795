#include <chrono>
#include <iostream>
#include <thread>
#include "Vector3.h"
#include "lodepng.h"
#include "timeutil.h"
#include "tinyexr.h"
#define THREAD_MULTIPLIER 1
void write_png(const std::vector<Vector3>& pixel_colors,
               const std::string& file_name, int width, int height);
void write_exr(const std::vector<Vector3>& hdr_image,
               const std::string& file_name, int width, int height);

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Please provide scene file as argument" << std::endl;
    return 1;
  }
  // Scene scene(argv[1]);
  std::cout << "Scene is parsed" << std::endl;
  const int thread_count =
      std::thread::hardware_concurrency() * THREAD_MULTIPLIER;

  /*const int camera_count = (int)scene.cameras.size();
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
    for (int j = 0; j < height; j++) {
      for (int i = 0; i < width; i++) {
        pixel_colors.push_back(pixels[j * width + i].get_color());
      }
    }

    std::string filename = camera.get_filename().substr(
        0, camera.get_filename().find_last_of("."));

    const Tonemapping_operator* tmo = camera.get_tmo();
    if (tmo) {
      write_exr(pixel_colors, filename, width, height);
      std::vector<Vector3> tonemapped_colors;
      tmo->apply_tmo(pixel_colors, tonemapped_colors);
      // TODO: parse gammacorrection style
      if (true) {
        constexpr float inverse = 1.0f / 2.4f;
        for (int j = 0; j < height; j++) {
          for (int i = 0; i < width; i++) {
            Vector3 old_rgb = tonemapped_colors[j * width + i];
            float r =
                clamp(
                    0.0f, 1.0f,
                    (1.055f * std::pow(clamp(0.0f, 1.0f, old_rgb.x), inverse) -
                     0.055f)) *
                255.0f;
            float g =
                clamp(
                    0.0f, 1.0f,
                    (1.055f * std::pow(clamp(0.0f, 1.0f, old_rgb.y), inverse) -
                     0.055f)) *
                255.0f;
            float b =
                clamp(
                    0.0f, 1.0f,
                    (1.055f * std::pow(clamp(0.0f, 1.0f, old_rgb.z), inverse) -
                     0.055f)) *
                255.0f;
            tonemapped_colors[j * width + i] = Vector3(r, g, b);
          }
        }
      }
      write_png(tonemapped_colors, filename, width, height);
    } else {
      write_png(pixel_colors, filename, width, height);
    }

    std::cout << filename << "(" << width << "x" << height << ") is saved in: ";
    print_time_diff(std::cout, start, end);
    std::cout << std::endl;
  }*/
  return 0;
}

void write_png(const std::vector<Vector3>& pixel_colors,
               const std::string& file_name, int width, int height) {
  unsigned char* image = new unsigned char[width * height * 4];
  int idx = 0;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      const Vector3& pixel = pixel_colors[j * width + i];
      image[idx++] = clamp(0.0f, 255.0f, pixel.x);
      image[idx++] = clamp(0.0f, 255.0f, pixel.y);
      image[idx++] = clamp(0.0f, 255.0f, pixel.z);
      image[idx++] = 255;
    }
  }
  std::string png_image_name = file_name + ".png";
  unsigned error =
      lodepng::encode(png_image_name.c_str(), image, width, height);
  // if there's an error, display it
  if (error) {
    std::cout << "encoder error " << error << ": " << lodepng_error_text(error)
              << std::endl;
  } else {
    printf("Saved png file. [ %s ] \n", png_image_name.c_str());
  }
  delete[] image;
}

void write_exr(const std::vector<Vector3>& hdr_image,
               const std::string& file_name, int width, int height) {
  EXRHeader header;
  InitEXRHeader(&header);

  EXRImage exr_image;
  InitEXRImage(&exr_image);

  exr_image.num_channels = 3;

  std::vector<float> images[3];
  images[0].resize(width * height);
  images[1].resize(width * height);
  images[2].resize(width * height);

  for (int i = 0; i < width * height; i++) {
    images[0][i] = hdr_image[i].x;
    images[1][i] = hdr_image[i].y;
    images[2][i] = hdr_image[i].z;
  }

  float* image_ptr[3];
  image_ptr[0] = &(images[2].at(0));  // B
  image_ptr[1] = &(images[1].at(0));  // G
  image_ptr[2] = &(images[0].at(0));  // R

  exr_image.images = (unsigned char**)image_ptr;
  exr_image.width = width;
  exr_image.height = height;

  header.num_channels = 3;
  header.channels =
      (EXRChannelInfo*)malloc(sizeof(EXRChannelInfo) * header.num_channels);
  // Must be (A)BGR order, since most of EXR viewers expect this
  // channel order.
  strncpy(header.channels[0].name, "B", 255);
  header.channels[0].name[strlen("B")] = '\0';
  strncpy(header.channels[1].name, "G", 255);
  header.channels[1].name[strlen("G")] = '\0';
  strncpy(header.channels[2].name, "R", 255);
  header.channels[2].name[strlen("R")] = '\0';

  header.pixel_types = (int*)malloc(sizeof(int) * header.num_channels);
  header.requested_pixel_types =
      (int*)malloc(sizeof(int) * header.num_channels);
  for (int i = 0; i < header.num_channels; i++) {
    header.pixel_types[i] =
        TINYEXR_PIXELTYPE_FLOAT;  // pixel type of input image
    header.requested_pixel_types[i] =
        TINYEXR_PIXELTYPE_HALF;  // pixel type of output image to
                                 // be stored in .EXR
  }

  std::string exr_image_name = file_name + ".exr";

  const char* err;
  int ret =
      SaveEXRImageToFile(&exr_image, &header, exr_image_name.c_str(), &err);
  if (ret != TINYEXR_SUCCESS) {
    fprintf(stderr, "Save EXR err: %s\n", err);
  }

  printf("Saved exr file. [ %s ] \n", exr_image_name.c_str());

  free(header.channels);
  free(header.pixel_types);
  free(header.requested_pixel_types);
}
