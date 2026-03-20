#
# Copyright (C) Rodolfo Herrera Hernandez. All rights reserved.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(opendxa VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release MinSizeRel RelWithDebInfo)

find_package(coretoolkit REQUIRED)
find_package(structure-identification REQUIRED)
find_package(common-neighbor-analysis REQUIRED)
find_package(polyhedral-template-matching REQUIRED)
find_package(TBB REQUIRED)
find_package(Boost REQUIRED)
find_package(spdlog REQUIRED)
find_package(nlohmann_json REQUIRED)

set(OPENDXA_BOOST_TARGET "")
if(TARGET Boost::headers)
	set(OPENDXA_BOOST_TARGET Boost::headers)
elseif(TARGET Boost::boost)
	set(OPENDXA_BOOST_TARGET Boost::boost)
elseif(TARGET boost::headers)
	set(OPENDXA_BOOST_TARGET boost::headers)
elseif(TARGET boost::boost)
	set(OPENDXA_BOOST_TARGET boost::boost)
else()
	message(FATAL_ERROR "OpenDXA requires a Boost headers target")
endif()

file(GLOB_RECURSE LIB_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/src/*.cpp)
list(FILTER LIB_SOURCES EXCLUDE REGEX ".*main\\.cpp$")

# Create static library with Volt code
add_library(${PROJECT_NAME}_lib STATIC ${LIB_SOURCES})
set_target_properties(${PROJECT_NAME}_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)

get_target_property(STRUCTURE_IDENTIFICATION_INCLUDE_DIRS structure-identification::structure-identification INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(COMMON_NEIGHBOR_ANALYSIS_INCLUDE_DIRS common-neighbor-analysis::common-neighbor-analysis INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(POLYHEDRAL_TEMPLATE_MATCHING_INCLUDE_DIRS polyhedral-template-matching::polyhedral-template-matching INTERFACE_INCLUDE_DIRECTORIES)

set(LOCAL_STRUCTURE_IDENTIFICATION_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/../StructureIdentification/include")
set(LOCAL_COMMON_NEIGHBOR_ANALYSIS_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/../CommonNeighborAnalysis/include")
set(LOCAL_POLYHEDRAL_TEMPLATE_MATCHING_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/../PolyhedralTemplateMatching/include")

if(EXISTS "${LOCAL_STRUCTURE_IDENTIFICATION_INCLUDE_DIR}")
	list(PREPEND STRUCTURE_IDENTIFICATION_INCLUDE_DIRS "${LOCAL_STRUCTURE_IDENTIFICATION_INCLUDE_DIR}")
endif()
if(EXISTS "${LOCAL_COMMON_NEIGHBOR_ANALYSIS_INCLUDE_DIR}")
	list(PREPEND COMMON_NEIGHBOR_ANALYSIS_INCLUDE_DIRS "${LOCAL_COMMON_NEIGHBOR_ANALYSIS_INCLUDE_DIR}")
endif()
if(EXISTS "${LOCAL_POLYHEDRAL_TEMPLATE_MATCHING_INCLUDE_DIR}")
	list(PREPEND POLYHEDRAL_TEMPLATE_MATCHING_INCLUDE_DIRS "${LOCAL_POLYHEDRAL_TEMPLATE_MATCHING_INCLUDE_DIR}")
endif()

target_include_directories(${PROJECT_NAME}_lib BEFORE PUBLIC
	${STRUCTURE_IDENTIFICATION_INCLUDE_DIRS}
	${COMMON_NEIGHBOR_ANALYSIS_INCLUDE_DIRS}
	${POLYHEDRAL_TEMPLATE_MATCHING_INCLUDE_DIRS}
	${CMAKE_SOURCE_DIR}/include
)
target_precompile_headers(${PROJECT_NAME}_lib PRIVATE <volt/core/volt.h>)
target_link_libraries(${PROJECT_NAME}_lib PUBLIC
	${OPENDXA_BOOST_TARGET}
	TBB::tbb
	structure-identification::structure-identification
	common-neighbor-analysis::common-neighbor-analysis
	polyhedral-template-matching::polyhedral-template-matching
	coretoolkit::coretoolkit
	spdlog::spdlog
	nlohmann_json::nlohmann_json
)

add_executable(${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/src/main.cpp)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_precompile_headers(${PROJECT_NAME} PRIVATE <volt/core/volt.h>)
target_link_libraries(${PROJECT_NAME} PRIVATE
	${PROJECT_NAME}_lib
	coretoolkit::coretoolkit
	TBB::tbb
)
install(TARGETS ${PROJECT_NAME} DESTINATION bin)

install(TARGETS ${PROJECT_NAME}_lib DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/ DESTINATION include)
