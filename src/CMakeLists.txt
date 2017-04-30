cmake_minimum_required(VERSION 3.1)

# get all cpp abs file names recursively starting from current directory
file(GLOB_RECURSE src_cpps *.cpp)
message(STATUS "Available src cpp files are: ${src_cpps}")

add_library(poa_alignment_lib STATIC ${src_cpps})
