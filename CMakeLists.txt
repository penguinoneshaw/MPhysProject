cmake_minimum_required(VERSION 3.12)
project(ProjectModelling VERSION 0.1.0)

add_compile_options(-Wextra -Wnon-virtual-dtor -pedantic)

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED COMPONENTS iostreams system filesystem)
find_package(Minuit2)
find_library(NCDFC NAMES netcdf_c++4 netcdf-cxx4)
find_library(NCDF NAMES netcdf libnetcdf)
include(CTest)
enable_testing()

add_executable(ProjectModelling main.cpp)
add_subdirectory(fitting)
add_subdirectory(grapher)

if (NOT NCDFC)
    add_subdirectory(netcdf-cxx4)
    set (NCDFC netcdf-cxx4)
endif()


set_property(TARGET ProjectModelling PROPERTY CXX_STANDARD 17)
set_property(TARGET ProjectModelling PROPERTY CXX_STANDARD_REQUIRED ON)
target_compile_features(ProjectModelling PRIVATE cxx_variable_templates cxx_lambdas cxx_lambda_init_captures )

if(HOME)
target_include_directories(ProjectModelling PRIVATE "${HOME}/include")
target_include_directories(ProjectFit PRIVATE "${HOME}/include")
endif(HOME)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

target_link_libraries(ProjectModelling PUBLIC ProjectFit)
target_link_libraries(ProjectModelling PRIVATE OpenMP::OpenMP_CXX ${Boost_LIBRARIES} Minuit2::Minuit2)

target_link_libraries(ProjectModelling PRIVATE ${NCDF})
target_link_libraries(ProjectModelling PRIVATE ${NCDFC})


