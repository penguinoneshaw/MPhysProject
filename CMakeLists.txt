cmake_minimum_required(VERSION 3.4.0)
project(ProjectModelling VERSION 0.1.0)

add_compile_options(-Wextra -Wnon-virtual-dtor -pedantic)

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED COMPONENTS iostreams system filesystem)
find_package(Minuit2)

include(CTest)
enable_testing()

add_executable(ProjectModelling main.cpp)
add_subdirectory(fitting)
add_subdirectory(grapher)


set_property(TARGET ProjectModelling PROPERTY CXX_STANDARD 14)
set_property(TARGET ProjectModelling PROPERTY CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(HOME)
target_include_directories(ProjectModelling PRIVATE "${HOME}/include")
target_include_directories(ProjectFit PRIVATE "${HOME}/include")
endif(HOME)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

target_link_libraries(ProjectModelling PUBLIC ProjectFit)
target_link_libraries(ProjectModelling PRIVATE netcdf-cxx4 netcdf OpenMP::OpenMP_CXX ${Boost_LIBRARIES} Minuit2::Minuit2)
