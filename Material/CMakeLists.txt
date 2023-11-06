cmake_minimum_required(VERSION 3.13.0)

project(Material)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../Libs)

include(Material1D/CMakeLists.txt)
include(Material2D/CMakeLists.txt)
include(Material3D/CMakeLists.txt)
include(MaterialOS/CMakeLists.txt)
include(Special/CMakeLists.txt)

add_library(${PROJECT_NAME} STATIC
        ExternalMaterial.cpp
        Material.cpp
        MaterialTester.cpp
        MaterialParser.cpp
        MaterialTemplate.cpp
        ${M1D}
        ${M2D}
        ${M3D}
        ${MOS}
        ${MS}
)
