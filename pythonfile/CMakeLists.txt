cmake_minimum_required(VERSION 3.20)
project(ProjectModellingPython VERSION 0.1.1 LANGUAGES CXX)

add_compile_options(-Wextra -Wnon-virtual-dtor -pedantic)

add_library(ProjectPython SHARED python_lib.cpp)

target_compile_features(ProjectPython PRIVATE cxx_std_17)
target_include_directories(ProjectPython PRIVATE ".." "../**")
