cmake_minimum_required(VERSION 3.13.0)

project(mumps C Fortran)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

include_directories(.)

add_compile_definitions(Add_)
add_compile_definitions(pord)
add_compile_definitions(metis)
add_compile_definitions(MUMPS_ARITH=MUMPS_ARITH_d)

include(../../Driver.cmake)

if ((CMAKE_C_COMPILER_ID MATCHES "GNU") OR (CMAKE_C_COMPILER_ID MATCHES "Clang")) # GNU GCC COMPILER
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
else () # MSVC COMPILER
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W0")
endif ()

file(GLOB mumps_f "*.F")
file(GLOB mumps_c "*.c")

add_library(mumps_f_mod OBJECT lr_common.F dlr_type.F dmumps_struc_def.F dlr_stats.F front_data_mgt_m.F dmumps_lr_data_m.F dlr_core.F dmumps_comm_buffer.F)
add_library(mumps_f_lib OBJECT ${mumps_f})
add_library(mumps_c_lib OBJECT ${mumps_c})

add_dependencies(mumps_f_lib mumps_f_mod)
add_dependencies(mumps_f_lib metis)
add_dependencies(mumps_c_lib metis)

add_library(${PROJECT_NAME} ${LIBRARY_TYPE} $<TARGET_OBJECTS:mumps_f_lib> $<TARGET_OBJECTS:mumps_c_lib>)

target_link_libraries(${PROJECT_NAME} metis)

message(STATUS "MUMPS Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
message(STATUS "MUMPS C_FLAGS: ${CMAKE_C_FLAGS}")