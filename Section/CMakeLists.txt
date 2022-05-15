cmake_minimum_required(VERSION 3.13.0)

project(Section)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../Libs)

include(Section1D/CMakeLists.txt)
include(Section2D/CMakeLists.txt)
include(Section3D/CMakeLists.txt)
include(SectionNM/CMakeLists.txt)
include(SectionShell/CMakeLists.txt)

add_library(${PROJECT_NAME} STATIC
        Section.cpp
        SectionParser.cpp
        ${S1D} ${S2D} ${S3D} ${SNM} ${SS}
        )