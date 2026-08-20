from conan import ConanFile
from conan.tools.cmake import CMake, cmake_layout


class ProjectModellingConan(ConanFile):
    required_conan_version = ">=2.28"
    name = "ProjectModelling"
    version = "1.0.0"
    description = "Physics project for modelling with CMake and Conan"
    author = "Your Name"
    settings = "os", "compiler", "build_type", "arch"
    options = {
        "shared": [True, False],
        "fPIC": [True, False],
        "simd": ["sse", "sse2", "avx", "avx2", "False"],
        "threads": [True, False],
    }
    default_options = {"shared": False, "fPIC": True, "simd": "avx2", "threads": True}
    requires = ("fftw/3.3.10", "netcdf/4.8.1")
    generators = "CMakeDeps", "CMakeToolchain"
    exports_sources = (
        "CMakeLists.txt",
        "src/*",
        "fitting/*",
        "grapher/*",
        "pythonfile/*",
        "*.hpp",
        "*.cpp",
    )

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def configure(self):
        # FFTW options
        self.options["fftw"].shared = self.options.shared
        self.options["fftw"].threads = self.options.threads
        self.options["fftw"].simd = self.options.simd
        self.options["fftw"].precision_longdouble = False
        

    def layout(self):
        cmake_layout(self)

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()
