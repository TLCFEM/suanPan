cmake_minimum_required(VERSION 3.13.0)

project(spmm Fortran)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

file(GLOB SPMM_SRC "*.f")

include(../../Driver.cmake)

add_library(spmm_obj OBJECT ${SPMM_SRC})

add_library(${PROJECT_NAME} ${LIBRARY_TYPE} $<TARGET_OBJECTS:spmm_obj>)

message(STATUS "SPMM Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")