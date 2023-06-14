cmake_minimum_required(VERSION 3.0)

project(lis C)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

set(CMAKE_C_STANDARD 11)

file(GLOB_RECURSE lis_src **/*.c)

include_directories(include)

add_compile_definitions(USE_VEC_COMP)

if (MSVC)
    add_compile_options(_CRT_SECURE_NO_WARNINGS)
endif ()

add_library(lis STATIC ${lis_src})

if (USE_MPI)
    add_library(lis-mpi STATIC ${lis_src})
    target_compile_definitions(lis-mpi PRIVATE USE_MPI)
endif ()
