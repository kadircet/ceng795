cmake_minimum_required (VERSION 3.11)

project (raytracer)
# The version number
set(raytracer_VERSION_MAJOR 3)
set(raytracer_VERSION_MINOR 0)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories("include")
include_directories("include/externals")
file(GLOB_RECURSE SOURCES "src/*.cpp")
file(GLOB_RECURSE HEADERS "include/*.h")

add_executable(raytracer ${SOURCES} ${HEADERS})
