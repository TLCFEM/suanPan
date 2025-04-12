cmake_minimum_required(VERSION 3.18)

project(superlu C)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

file(GLOB SUPERLU_SRC "*.c")

add_compile_definitions(NDEBUG)
add_compile_definitions(PRNTlevel=0)
add_compile_definitions(Add_)

if (CMAKE_C_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W0")
else ()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
endif ()

add_library(${PROJECT_NAME} STATIC ${SUPERLU_SRC})

message(STATUS "SUPERLU C_FLAGS: ${CMAKE_C_FLAGS}")