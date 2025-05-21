cmake_minimum_required(VERSION 3.18)

project(MaterialExample CXX)

set(CMAKE_CXX_STANDARD 20)

add_compile_definitions(SUANPAN_DLL)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

link_directories(${LIBRARY_OUTPUT_PATH}/../)

if (CMAKE_SYSTEM_NAME MATCHES "Darwin") # MAC PLATFORM
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,-U,_SUANPAN_PRINT")
endif ()

add_library(${PROJECT_NAME} SHARED
        ElasticExternal.cpp
        MaterialExample.cpp
        ../../Domain/Tag.cpp
        ../../Material/Material.cpp
        ../../Include/fmt/src/format.cc
)

target_link_libraries(${PROJECT_NAME} suanPan)

install(TARGETS ${PROJECT_NAME})

message(STATUS "${PROJECT_NAME} CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
