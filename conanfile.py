from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, CMakeDeps, cmake_layout


class OpenDXAConan(ConanFile):
    name = "opendxa"
    version = "1.0.0"
    package_type = "application"
    license = "MIT"
    settings = "os", "arch", "compiler", "build_type"
    default_options = {
        "hwloc/*:shared": True,
    }
    requires = (
        "boost/1.88.0",
        "onetbb/2021.12.0",
        "coretoolkit/1.0.0",
        "structure-identification/1.0.0",
        "spdlog/1.14.1",
        "nlohmann_json/3.11.3",
    )
    exports_sources = "CMakeLists.txt", "include/*", "src/*"

    def layout(self):
        cmake_layout(self)

    def generate(self):
        CMakeToolchain(self).generate()
        CMakeDeps(self).generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()
