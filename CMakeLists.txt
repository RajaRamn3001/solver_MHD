cmake_minimum_required(VERSION 3.5)

project(myProject)

enable_language(Fortran)

message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/run)

add_subdirectory(source)
