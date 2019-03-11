# "Singing in a Warming Ocean" - SOFAR Channel Movements

## Building
The project is set up to be built with cmake 3.12 or greater, and requires a C++17 compatible compiler (which may require running `module add devtoolset/7` before configuring). It is designed to be run on the Met Office EN4 profile data dataset, but will work on any dataset which uses 
- LAT
- LONG
- DEPH_CORRECTED
- POTM_CORRECTED
- PSAL_CORRECTED

as its indexed variables. These can relatively easily be changed in the `main.cpp` file.

```zsh
mkdir build
cd build
cmake3 .. -D CMAKE_BUILD_PARALLEL_LEVEL=4
cmake3 --build . --config Release -j 4
cd ..
build/ProjectModelling <directory name containing data>
```

## Other Build Requirements
 - the `netcdf` C library
 - Boost (with the up-to-date filesystems module)
 - 'Fastest Fourier Transform in the West' (fftw3)