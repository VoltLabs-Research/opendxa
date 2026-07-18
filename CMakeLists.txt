cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(opendxa VERSION 1.0.0 LANGUAGES CXX)

find_package(coretoolkit REQUIRED)
find_package(structure-identification REQUIRED)
find_package(yaml-cpp REQUIRED)  # transitive dep of structure-identification (static .a), must be linked into the exe

volt_add_plugin(opendxa
    MAIN_SOURCE src/plugin.cpp
    DEPENDENCIES
        structure-identification::structure-identification
        yaml-cpp::yaml-cpp
    INSTALL_DATA lattices
)
