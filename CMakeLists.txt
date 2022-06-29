cmake_minimum_required(VERSION 3.12)
project(ProjectModelling VERSION 1.0.0)

add_compile_options(-Wextra -Wnon-virtual-dtor -pedantic -fconcepts)

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED COMPONENTS iostreams system filesystem)

add_subdirectory(Minuit2)

find_library(NCDFC NAMES netcdf_c++4 netcdf-cxx4)
find_library(NCDF NAMES netcdf libnetcdf)

add_executable(ProjectModelling main.cpp)
add_subdirectory(fitting)
add_subdirectory(grapher)

if (NOT NCDFC)
    add_subdirectory(netcdf-cxx4)
    set (NCDFC netcdf-cxx4)
endif()


set_target_properties(ProjectModelling PROPERTIES CXX_STANDARD 20 
    CXX_STANDARD_REQUIRED TRUE 
    CXX_EXTENSIONS TRUE
    )

target_compile_features(ProjectModelling PRIVATE cxx_variable_templates cxx_lambdas cxx_lambda_init_captures cxx_digit_separators)
target_compile_options(ProjectModelling PRIVATE -fconcepts)

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


