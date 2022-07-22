cmake_minimum_required(VERSION 3.13.0)

project(metis C)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

if (CMAKE_C_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W0")
else ()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
endif ()

include(gklib/CMakeLists.txt)
include_directories(gklib)

add_subdirectory(src)
#add_subdirectory(programs)
