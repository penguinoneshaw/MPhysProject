add_library(ProjectFit STATIC projectfit.cpp)

set_property(TARGET ProjectFit PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET ProjectFit PROPERTY CXX_STANDARD 17)

target_link_libraries(ProjectFit PUBLIC fftw3f fftw3f_threads fftw3 fftw3_threads)
target_include_directories(ProjectFit PRIVATE "../")

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED COMPONENTS iostreams system filesystem)
target_link_libraries(ProjectFit PUBLIC OpenMP::OpenMP_CXX ${Boost_LIBRARIES} Minuit2::Minuit2)
