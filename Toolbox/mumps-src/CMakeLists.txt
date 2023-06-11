cmake_minimum_required(VERSION 3.13.0)

project(mumps C Fortran)

# ifx does not compile properly due to potential memory allocation issues
if (CMAKE_Fortran_COMPILER MATCHES "ifx")
    set(CMAKE_Fortran_COMPILER "ifort")
endif ()

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

include_directories(.)

add_compile_definitions(Add_)
add_compile_definitions(pord)
add_compile_definitions(metis)

if (USE_MKL)
    add_compile_definitions(GEMMT_AVAILABLE)
endif ()

include(../../Driver.cmake)

if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W0")
else ()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
endif ()

file(GLOB mumps_f "*.F")
file(GLOB mumps_c "*.c")

set(mumps_mod
        lr_common.F
        dlr_type.F
        dmumps_struc_def.F
        dlr_stats.F
        slr_type.F
        smumps_struc_def.F
        slr_stats.F
        front_data_mgt_m.F
        dmumps_lr_data_m.F
        dlr_core.F
        dmumps_comm_buffer.F
        smumps_lr_data_m.F
        slr_core.F
        smumps_comm_buffer.F
        )
set(mumps_mpi_c elapse.c mpic.c)
set(mumps_mpi_f mpi.F)

foreach (mod_src ${mumps_mod})
    foreach (f_src ${mumps_f})
        if (f_src MATCHES ${mod_src})
            list(REMOVE_ITEM mumps_f ${f_src})
            break()
        endif ()
    endforeach ()
endforeach ()

foreach (mod_src ${mumps_mpi_f})
    foreach (f_src ${mumps_f})
        if (f_src MATCHES ${mod_src})
            list(REMOVE_ITEM mumps_f ${f_src})
            break()
        endif ()
    endforeach ()
endforeach ()

foreach (mod_src ${mumps_mpi_c})
    foreach (c_src ${mumps_c})
        if (c_src MATCHES ${mod_src})
            list(REMOVE_ITEM mumps_c ${c_src})
            break()
        endif ()
    endforeach ()
endforeach ()

set_source_files_properties(dmumps_c.c PROPERTIES COMPILE_DEFINITIONS "MUMPS_ARITH=MUMPS_ARITH_d")
set_source_files_properties(smumps_c.c PROPERTIES COMPILE_DEFINITIONS "MUMPS_ARITH=MUMPS_ARITH_s")

add_library(mumps_mod_lib OBJECT ${mumps_mod})
add_library(mumps_f_lib OBJECT ${mumps_f})
add_library(mumps_c_lib OBJECT ${mumps_c})

add_dependencies(mumps_f_lib mumps_mod_lib)
add_dependencies(mumps_f_lib metis)
add_dependencies(mumps_c_lib metis)

if (USE_MPI)
    add_library(${PROJECT_NAME} ${LIBRARY_TYPE}
            $<TARGET_OBJECTS:mumps_mod_lib>
            $<TARGET_OBJECTS:mumps_f_lib>
            $<TARGET_OBJECTS:mumps_c_lib>)
else ()
    add_library(mumps_mpi_f_lib OBJECT ${mumps_mpi_f})
    add_library(mumps_mpi_c_lib OBJECT ${mumps_mpi_c})
    add_library(${PROJECT_NAME} ${LIBRARY_TYPE}
            $<TARGET_OBJECTS:mumps_mod_lib>
            $<TARGET_OBJECTS:mumps_mpi_f_lib>
            $<TARGET_OBJECTS:mumps_mpi_c_lib>
            $<TARGET_OBJECTS:mumps_f_lib>
            $<TARGET_OBJECTS:mumps_c_lib>)
endif ()

target_link_libraries(${PROJECT_NAME} metis)

message(STATUS "MUMPS Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
message(STATUS "MUMPS C_FLAGS: ${CMAKE_C_FLAGS}")
