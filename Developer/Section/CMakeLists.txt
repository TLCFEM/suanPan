cmake_minimum_required(VERSION 3.0)

project(SectionExample CXX)

set(CMAKE_CXX_STANDARD 17)

add_definitions(-DSUANPAN_DLL)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

link_directories(${LIBRARY_OUTPUT_PATH}/../)

if(CMAKE_SYSTEM_NAME MATCHES "Darwin") # MAC PLATFORM
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,-U,_SUANPAN_PRINT")
endif()

if(COMPILER_IDENTIFIER MATCHES "vs")
    link_libraries(libopenblas)
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    link_libraries(openblas)
endif()

link_libraries(suanPan)

add_library(${PROJECT_NAME} SHARED
    SectionExample.cpp
    ../../Domain/Tag.cpp
    ../../Material/Material.cpp
    ../../Section/Section.cpp
    ../../Toolbox/debug.cpp
    )

message("SectionExample CXX_FLAGS: ${CMAKE_CXX_FLAGS}")