cmake_minimum_required(VERSION 3.18)

project(lis C)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

set(CMAKE_C_STANDARD 11)

file(GLOB_RECURSE lis_src **/*.c)

include_directories(include)

add_compile_definitions(USE_VEC_COMP)

if (MSVC)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
endif ()

if (CMAKE_C_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W0")
else ()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
endif ()

add_library(lis STATIC ${lis_src})

message(STATUS "LIS C_FLAGS: ${CMAKE_C_FLAGS}")
