cmake_minimum_required(VERSION 3.18)

project(amd Fortran)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

file(GLOB AMD_SRC "*.f")

include(../../Driver.cmake)

add_library(amd_obj OBJECT ${AMD_SRC})

add_library(${PROJECT_NAME} ${LIBRARY_TYPE} $<TARGET_OBJECTS:amd_obj>)

message(STATUS "AMD Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")