cmake_minimum_required(VERSION 3.4.0)
project(ProjectModelling VERSION 0.1.0)

add_compile_options(-Wextra -Wnon-virtual-dtor -pedantic)

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED COMPONENTS iostreams system filesystem)
find_package(Minuit2)
find_library(NCDF netcdf_c++4)
find_library(NCDFC netcdf)
include(CTest)
enable_testing()

add_executable(ProjectModelling main.cpp)
add_subdirectory(fitting)
add_subdirectory(grapher)


set_property(TARGET ProjectModelling PROPERTY CXX_STANDARD 14)
set_property(TARGET ProjectModelling PROPERTY CXX_STANDARD_REQUIRED ON)
target_compile_features(ProjectModelling PRIVATE cxx_variable_templates)

if(HOME)
target_include_directories(ProjectModelling PRIVATE "${HOME}/include")
target_include_directories(ProjectFit PRIVATE "${HOME}/include")
endif(HOME)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

target_link_libraries(ProjectModelling PUBLIC ProjectFit)
target_link_libraries(ProjectModelling PRIVATE ${NCDF} ${NCDFC} OpenMP::OpenMP_CXX ${Boost_LIBRARIES} Minuit2::Minuit2)
