# Project Dependencies and Resolution Strategy

## Required Dependencies

### 1. FFTW 3.3.10 (Fast Fourier Transform Library)
- **Source**: Conan (ConanCenter)
- **Status**: Required, managed via Conan
- **CMake Target**: `FFTW::fftw3f`, `FFTW::fftw3f_threads`, `FFTW::fftw3`, `FFTW::fftw3_threads`
- **Used by**: `ProjectFit` library

### 2. NetCDF-C 4.9.2 (Data I/O Library)
- **Source**: Conan (ConanCenter)
- **Status**: Required, managed via Conan
- **CMake Target**: `netCDF::netcdf`
- **Used by**: Main executable (`ProjectModelling`)

### 3. Minuit2 (Minimization Library)
- **Status**: Required
- **Resolution Strategy** (in priority order):
  1. **Check for ROOT**: If ROOT 6.x is found on the system via `find_package(ROOT)`, use `ROOT::Minuit2`
  2. **Build from submodule**: If ROOT is not found, build Minuit2 from the submodule at `/Minuit2/`
- **Used by**: `ProjectFit` library for curve fitting operations

## Build Process

### Step 1: Conan Install
```bash
conan install . --build=missing
```
- Downloads and/or builds FFTW and NetCDF-C
- Generates CMake configuration files in `build/Release/generators/`

### Step 2: CMake Configure
```bash
cmake --preset conan-release
```
or
```bash
cd build/Release
cmake ../.. -DCMAKE_TOOLCHAIN_FILE=generators/conan_toolchain.cmake
```

**What happens during configure:**
- Finds FFTW and NetCDF-C via Conan-generated find modules
- Looks for ROOT on the system
  - If found: Uses ROOT's Minuit2
  - If not found: Configures Minuit2 submodule for build
- Configures all project targets (ProjectModelling, ProjectFit, Grapher)

### Step 3: CMake Build
```bash
cmake --build build/Release
```

**Build order** (automatic dependency resolution):
1. Minuit2 (if building from submodule)
2. FFTW and NetCDF-C (already built by Conan)
3. ProjectFit (depends on FFTW and Minuit2)
4. Grapher
5. ProjectModelling (main executable)

## Compiler Requirements

- **C++ Standard**: C++23 (enforced via `target_compile_features()`)
- **Compiler**: 
  - GCC 13+
  - Clang 17+
  - MSVC 2022+

## Troubleshooting

### If Conan fails to find FFTW or NetCDF-C
```bash
conan remote list  # Verify remotes are configured
conan search fftw/*  # Check if package exists
conan install . --build=missing -v  # Verbose mode for debugging
```

### If Minuit2 build fails
- Check submodule is initialized: `git submodule update --init --recursive`
- Check compiler version: `gcc --version` or `clang --version`
- Ensure C++23 support in compiler
- Build with verbose output: `cmake --build build/Release --verbose`

### If ROOT is found but you want to use the submodule
Edit the main `CMakeLists.txt` and remove or comment out:
```cmake
find_package(ROOT QUIET CONFIG)
```
