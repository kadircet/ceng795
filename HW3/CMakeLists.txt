cmake_minimum_required (VERSION 3.11)

project (raytracer)
# The version number
set(raytracer_VERSION_MAJOR 3)
set(raytracer_VERSION_MINOR 0)

include_directories(include)
file(GLOB SOURCES "src/*.cpp")

add_executable(raytracer ${SOURCES})