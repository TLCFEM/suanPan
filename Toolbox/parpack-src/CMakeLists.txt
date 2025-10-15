cmake_minimum_required(VERSION 3.18)

project(parpack Fortran)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

file(GLOB PARPACK_SRC "*.f")

add_library(parpack_obj OBJECT ${PARPACK_SRC})

add_library(${PROJECT_NAME} $<TARGET_OBJECTS:parpack_obj>)

target_link_libraries(${PROJECT_NAME} arpack MPI::MPI_Fortran)

message(STATUS "PARPACK Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")