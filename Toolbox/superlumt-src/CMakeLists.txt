cmake_minimum_required(VERSION 3.18)

project(superlumt C)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

include(../../Driver.cmake)

file(GLOB SUPERLUMT_SRC "*.c")

add_compile_definitions(NDEBUG)
add_compile_definitions(PRNTlevel=0)
add_compile_definitions(Add_)

if (CMAKE_C_COMPILER_ID MATCHES "MSVC")
    add_compile_definitions(__OPENMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W0 /openmp")
    set_source_files_properties(await.c PROPERTIES COMPILE_FLAGS "/Od")
else ()
    add_compile_definitions(__PTHREAD)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
    set_source_files_properties(await.c PROPERTIES COMPILE_FLAGS "-O0")
endif ()

add_library(${PROJECT_NAME} ${LIBRARY_TYPE} ${SUPERLUMT_SRC})

message(STATUS "SUPERLUMT C_FLAGS: ${CMAKE_C_FLAGS}")