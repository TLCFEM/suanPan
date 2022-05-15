cmake_minimum_required(VERSION 3.13.0)

project(metis C)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

include(gklib/CMakeLists.txt)
include_directories(gklib)

add_subdirectory(src)
#add_subdirectory(programs)
