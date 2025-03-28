cmake_minimum_required(VERSION 3.18)

project(metis C)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

if (USE_64BIT_INDEXING) 
    add_compile_definitions(IDXTYPEWIDTH=64)
    message(STATUS "Using 64-bit indexing for METIS")
endif()

if (CMAKE_C_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W0")
else ()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
endif ()

include(gklib/CMakeLists.txt)
include_directories(gklib)

add_subdirectory(src)
#add_subdirectory(programs)
