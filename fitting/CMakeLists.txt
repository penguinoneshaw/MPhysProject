cmake_minimum_required(VERSION 3.20)
project(ProjectFitting VERSION 1.0.0 LANGUAGES CXX)

# Required dependencies
find_package(FFTW3 REQUIRED CONFIG)
find_package(ROOT QUIET CONFIG)
if(NOT ROOT_FOUND)
  # Minuit2 should be available from parent project's subdirectory
endif()

add_library(ProjectFit STATIC
  projectfit.cpp
  ideal_fit_function.cpp
  quadratic_fit_function.cpp
)

target_compile_features(ProjectFit PRIVATE cxx_std_23)
target_include_directories(ProjectFit PRIVATE "../")

# Link FFTW from Conan (required)
target_link_libraries(ProjectFit PUBLIC fftw::fftw)

# Link Minuit2 (required - from ROOT or submodule)
if(ROOT_FOUND)
  target_link_libraries(ProjectFit PUBLIC ROOT::Minuit2)
else()
  target_link_libraries(ProjectFit PUBLIC Minuit2)
endif()
