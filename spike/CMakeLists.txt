cmake_minimum_required(VERSION 3.1)

# get all cpp abs file names recursively starting from current directory
file(GLOB_RECURSE spike_cpps *.cpp)
message(STATUS "Available spike cpp files are: ${spike_cpps}")

add_library(spike_lib STATIC ${spike_cpps})
