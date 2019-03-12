cmake_minimum_required(VERSION 3.12)
project(ProjectModellingPython VERSION 0.1.1)
add_compile_options(-Wextra -Wnon-virtual-dtor -pedantic)

add_library(ProjectPython SHARED python_lib.cpp)

target_include_directories(ProjectPython PRIVATE "..")
set_property(TARGET ProjectPython PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET ProjectPython PROPERTY CXX_STANDARD 17)