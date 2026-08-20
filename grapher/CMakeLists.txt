cmake_minimum_required(VERSION 3.20)
project(Grapher VERSION 1.0.0 LANGUAGES CXX)

# Create a header-only library interface
add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
