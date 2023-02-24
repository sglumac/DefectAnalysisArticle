cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 11)

set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

project(Defects)

add_subdirectory(src)
add_subdirectory(scripts)
add_subdirectory(tests)
