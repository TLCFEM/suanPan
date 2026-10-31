cmake_minimum_required(VERSION 3.18)

project(replxx CXX)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

add_library(${PROJECT_NAME} STATIC
        ConvertUTF.cpp
        conversion.cxx
        escape.cxx
        history.cxx
        prompt.cxx
        replxx.cxx
        replxx_impl.cxx
        terminal.cxx
        util.cxx
        wcwidth.cpp
        windows.cxx
)

if (WIN32)
    add_compile_definitions(REPLXX_STATIC)
endif ()

if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    string(REGEX REPLACE "/we(4244|4267|4305)" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0")
else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
endif ()

find_package(Threads REQUIRED)

link_libraries(Threads::Threads)

message(STATUS "REPLXX CXX_FLAGS: ${CMAKE_CXX_FLAGS}")