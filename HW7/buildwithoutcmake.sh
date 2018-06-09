g++ src/*.cpp -o raytracer $(pkg-config --cflags --libs opencv) -std=c++14 -I include/ -O3 
