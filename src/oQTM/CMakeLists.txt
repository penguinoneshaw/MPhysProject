cmake_minimum_required(VERSION 3.20)
project(oQTM LANGUAGES CXX)

add_library(oQTM INTERFACE)
target_include_directories(oQTM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(oQTM INTERFACE cxx_std_23)

find_package(Threads REQUIRED)
target_link_libraries(oQTM INTERFACE Threads::Threads)