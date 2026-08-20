cmake_minimum_required(VERSION 3.20)
project(ProjectModelling VERSION 1.0.0 LANGUAGES CXX)

# Compiler warnings
add_compile_options(-Wextra -Wall -Wnon-virtual-dtor -Wpedantic)

# Address sanitizer for debugging (optional, uncomment to enable)
# add_compile_options(-fsanitize=address)
# add_link_options(-fsanitize=address)

# Find Conan-managed dependencies (required)
find_package(FFTW3 REQUIRED CONFIG)
find_package(netCDF REQUIRED CONFIG)
find_package(Boost REQUIRED CONFIG)

# Find or build Minuit2 (required)
find_package(ROOT QUIET CONFIG)
if(NOT ROOT_FOUND)
  # Build Minuit2 from submodule if ROOT not found
  message(STATUS "ROOT not found, building Minuit2 from submodule")
  add_subdirectory(external/Minuit2)
else()
  message(STATUS "ROOT found at ${ROOT_DIR}")
endif()

add_library(gnuplot-iostream INTERFACE)
target_include_directories(gnuplot-iostream INTERFACE ${CMAKE_SOURCE_DIR}/external/gnuplot-iostream)

find_library(netCDFCXX REQUIRED NAMES netcdf_c++4 netcdf-c++4  netcdf-cxx4)


# Add executable
add_executable(ProjectModelling src/main.cpp)
target_compile_features(ProjectModelling PRIVATE cxx_std_23)

# Add subdirectories
add_subdirectory(src/grapher)
add_subdirectory(src/fitting)
add_subdirectory(src/oQTM)

# Link libraries
target_link_libraries(ProjectModelling PRIVATE
  ProjectFit
  Grapher
  oQTM
  netCDF::netcdf
  ${netCDFCXX}
)

# CPack configuration
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

enable_testing()
find_package(GTest REQUIRED)
add_subdirectory(tests)