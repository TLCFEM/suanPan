cmake_minimum_required(VERSION 3.18)

project(fext Fortran)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

add_library(${PROJECT_NAME}
        $<TARGET_OBJECTS:amd_obj>
        $<TARGET_OBJECTS:arpack_obj>
        $<TARGET_OBJECTS:feast_obj>
        $<TARGET_OBJECTS:mumps_mpi_f_lib>
        $<TARGET_OBJECTS:mumps_mpi_c_lib>
        $<TARGET_OBJECTS:mumps_f_lib>
        $<TARGET_OBJECTS:mumps_c_lib>
)

target_link_libraries(${PROJECT_NAME} metis)
