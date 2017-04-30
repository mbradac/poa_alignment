cmake_minimum_required(VERSION 3.1)

# set project name
project(poa_alignment)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++14")

# dir variables
set(src_dir ${CMAKE_SOURCE_DIR}/src)
set(libs_dir ${CMAKE_SOURCE_DIR}/libs)
set(tests_dir ${CMAKE_SOURCE_DIR}/tests)
set(spike_dir ${CMAKE_SOURCE_DIR}/spike)

# release flags
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
#debug flags
set(CMAKE_CXX_FLAGS_DEBUG "-g")

if ("${CMAKE_BUILD_TYPE}" STREQUAL "")
    set(CMAKE_BUILD_TYPE "debug")
endif()
message(STATUS "CMake build type: ${CMAKE_BUILD_TYPE}")

include_directories(${src_dir})
include_directories(${libs_dir})

add_subdirectory(${src_dir})
add_subdirectory(${spike_dir})
add_subdirectory(${tests_dir})
