cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(opendxa VERSION 1.0.0 LANGUAGES CXX)

find_package(coretoolkit REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(structure-identification REQUIRED)

file(GLOB PTM_MATH_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/dependencies/ptm_math/*.cpp)

volt_add_plugin(opendxa
    MAIN_SOURCE src/plugin.cpp
    VENDOR_SOURCES ${PTM_MATH_SOURCES}
    INCLUDE_DIRS dependencies/ptm_math
    DEPENDENCIES
        structure-identification::structure-identification
        yaml-cpp::yaml-cpp
    INSTALL_DATA lattices
)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/structures/
    DESTINATION share/volt/structures)
