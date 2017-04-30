cmake_minimum_required(VERSION 3.1)

project(tests)

enable_testing()

# get all cpp abs file names recursively starting from current directory
file(GLOB_RECURSE test_cpps *.cpp)
message(STATUS "Available tests cpp files are: ${test_cpps}")

# for each cpp file build binary and register test
foreach(test_cpp ${test_cpps})
    # get exec name (remove extension from the abs path)
    get_filename_component(exec_name ${test_cpp} NAME_WE)

    # set target name
    set(target_name ${exec_name})

    # build exec file
    add_executable(${target_name} ${test_cpp})
    set_property(TARGET ${target_name} PROPERTY CXX_STANDARD ${cxx_standard})

    # link libraries
    target_link_libraries(${target_name} poa_alignment_lib)
    target_link_libraries(${target_name} spike_lib)

    # register test
    add_test(${target_name} ${exec_name})
endforeach()

# copy test data
file(COPY ${CMAKE_SOURCE_DIR}/tests/data
     DESTINATION ${CMAKE_BINARY_DIR}/tests)
