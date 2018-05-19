cmake_minimum_required (VERSION 3.11)

project (raytracer)
# The version number
set(raytracer_VERSION_MAJOR 6)
set(raytracer_VERSION_MINOR 0)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories("include")
file(GLOB_RECURSE SOURCES "src/*.cpp")
file(GLOB_RECURSE HEADERS "include/*.h")
# on windows, download opencv from https://sourceforge.net/projects/opencvlibrary/files/opencv-win/3.4.1/opencv-3.4.1-vc14_vc15.exe/download
# after extraction, copy files under build folder to $ENV{ProgramFiles}/OpenCV, and after running this script
# console will tell you a path which you should add to Path variable
find_package( OpenCV REQUIRED )
add_executable(raytracer ${SOURCES} ${HEADERS})
target_link_libraries( raytracer ${OpenCV_LIBS} )
