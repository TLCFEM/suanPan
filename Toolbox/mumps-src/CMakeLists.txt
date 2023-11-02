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
        mumps_pivnul_mod.F
        ana_omp_m.F
        ana_blk_m.F
        ana_orderings_wrappers_m.F
        double_linked_list.F
        fac_asm_build_sort_index_ELT_m.F
        fac_asm_build_sort_index_m.F
        fac_descband_data_m.F
        fac_future_niv2_mod.F
        fac_ibct_data_m.F
        fac_maprow_data_m.F
        front_data_mgt_m.F
        lr_common.F
        mumps_comm_ibcast.F
        mumps_l0_omp_m.F
        mumps_memory_mod.F
        mumps_mpitoomp_m.F
        mumps_ooc_common.F
        mumps_static_mapping.F
        omp_tps_common_m.F

        dana_aux.F
        dana_aux_par.F
        dana_lr.F
        dfac_asm_master_ELT_m.F
        dfac_asm_master_m.F
        dfac_front_aux.F
        dfac_front_LU_type1.F
        dfac_front_LU_type2.F
        dfac_front_LDLT_type1.F
        dfac_front_LDLT_type2.F
        dfac_front_type2_aux.F
        dfac_sispointers_m.F
        dfac_lr.F
        dfac_mem_dynamic.F
        dfac_omp_m.F
        dfac_par_m.F
        dfac_sol_l0omp_m.F
        dlr_core.F
        dlr_stats.F
        dlr_type.F
        dmumps_comm_buffer.F
        dmumps_config_file.F
        dmumps_load.F
        dmumps_lr_data_m.F
        dmumps_ooc_buffer.F
        dmumps_ooc.F
        dmumps_sol_es.F
        dmumps_save_restore.F
        dmumps_save_restore_files.F
        dmumps_struc_def.F
        domp_tps_m.F
        dsol_lr.F
        dsol_omp_m.F
        dstatic_ptr_m.F
        dmumps_mpi3_mod.F

        sana_aux.F
        sana_aux_par.F
        sana_lr.F
        sfac_asm_master_ELT_m.F
        sfac_asm_master_m.F
        sfac_front_aux.F
        sfac_front_LU_type1.F
        sfac_front_LU_type2.F
        sfac_front_LDLT_type1.F
        sfac_front_LDLT_type2.F
        sfac_front_type2_aux.F
        sfac_sispointers_m.F
        sfac_lr.F
        sfac_mem_dynamic.F
        sfac_omp_m.F
        sfac_par_m.F
        sfac_sol_l0omp_m.F
        slr_core.F
        slr_stats.F
        slr_type.F
        smumps_comm_buffer.F
        smumps_config_file.F
        smumps_load.F
        smumps_lr_data_m.F
        smumps_ooc_buffer.F
        smumps_ooc.F
        smumps_sol_es.F
        smumps_save_restore.F
        smumps_save_restore_files.F
        smumps_struc_def.F
        somp_tps_m.F
        ssol_lr.F
        ssol_omp_m.F
        sstatic_ptr_m.F
        smumps_mpi3_mod.F
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
