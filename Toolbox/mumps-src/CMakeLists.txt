cmake_minimum_required(VERSION 3.0)

project(mumps C Fortran)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

include_directories(.)

add_definitions(-DAdd_ -Dpord -DMUMPS_ARITH=MUMPS_ARITH_d)

if(CMAKE_C_COMPILER_ID MATCHES "GNU") # GNU GCC COMPILER
    link_libraries(pthread gfortran quadmath openblas)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
else() # MSVC COMPILER
    link_libraries(libopenblas)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W0")
endif()

file(GLOB mumps_f "*.F")
file(GLOB mumps_c "*.c")

add_library(mumps_f_mod OBJECT dlr_type.F dmumps_comm_buffer.F)
add_library(mumps_f_lib OBJECT ${mumps_f})
add_library(mumps_c_lib OBJECT ${mumps_c})

add_dependencies(mumps_f_lib mumps_f_mod)

if(BUILD_SHARED)
    add_library(${PROJECT_NAME} SHARED $<TARGET_OBJECTS:mumps_f_lib> $<TARGET_OBJECTS:mumps_c_lib>)
else()
    add_library(${PROJECT_NAME} STATIC $<TARGET_OBJECTS:mumps_f_lib> $<TARGET_OBJECTS:mumps_c_lib>)
endif()

message("MUMPS Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
message("MUMPS C_FLAGS: ${CMAKE_C_FLAGS}")