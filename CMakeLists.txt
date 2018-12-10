cmake_minimum_required(VERSION 3.4.0)
project(ProjectModelling VERSION 0.1.0)

add_compile_options(-Wextra -Wshadow -Wnon-virtual-dtor -pedantic)

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED COMPONENTS iostreams system filesystem)

find_library(NETCDFCPP netcdf-cxx4)

include(CTest)
enable_testing()

add_executable(ProjectModelling main.cpp)
target_link_libraries(ProjectModelling PRIVATE OpenMP::OpenMP_CXX)
add_subdirectory(grapher)
target_link_libraries(ProjectModelling PUBLIC ${Boost_LIBRARIES})


set_property(TARGET ProjectModelling PROPERTY CXX_STANDARD 14)
set_property(TARGET ProjectModelling PROPERTY CXX_STANDARD_REQUIRED ON)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

set_property(TARGET ProjectModelling PROPERTY CXX_STANDARD 17)
target_link_libraries( ${PROJECT_NAME} PUBLIC ${NETCDFCPP} fftw3f fftw3f_threads)



