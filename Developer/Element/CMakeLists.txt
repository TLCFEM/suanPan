cmake_minimum_required(VERSION 3.18)

project(ElementExample CXX)

set(CMAKE_CXX_STANDARD 20)

add_compile_definitions(SUANPAN_DLL)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

link_directories(${LIBRARY_OUTPUT_PATH}/../)

if (CMAKE_SYSTEM_NAME MATCHES "Darwin") # MAC PLATFORM
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,-U,_SUANPAN_PRINT")
endif ()

link_libraries(suanPan)

add_library(${PROJECT_NAME} SHARED
        ElementExample.cpp
        ../../Domain/Group/Group.cpp
        ../../Domain/Node.cpp
        ../../Domain/Tag.cpp
        ../../Element/Element.cpp
        ../../Material/Material.cpp
        ../../Section/Section.cpp
        ../../Include/fmt/src/format.cc
)

add_dependencies(ElementExample suanPan)

message(STATUS "ElementExample CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
