cmake_minimum_required(VERSION 3.18)

project(feast Fortran)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

file(GLOB FEAST_SRC *.f90)

set(SPIKE_SRC
        spike-src/cpivul.f
        spike-src/cspike_smp.f90
        spike-src/dpivul.f
        spike-src/dspike_smp.f90
        spike-src/dzlbprim.f90
        spike-src/f90_functions_wrapper.f90
        spike-src/sclbprim.f90
        spike-src/spike_smp_utilities.f90
        spike-src/spivul.f
        spike-src/sspike_smp.f90
        spike-src/zpivul.f
        spike-src/zspike_smp.f90
)

if (SP_ENABLE_MKL)
    add_compile_definitions(MKL)
endif ()

if (NOT COMPILER_IDENTIFIER MATCHES "vs") # GNU GCC COMPILER
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -ffixed-line-length-none")
endif ()

add_library(feast_obj OBJECT ${FEAST_SRC} ${SPIKE_SRC})

add_library(${PROJECT_NAME} $<TARGET_OBJECTS:feast_obj>)

if (SP_ENABLE_MPI)
    add_library(feast_mpi_obj OBJECT ${FEAST_SRC} ${SPIKE_SRC})
    target_compile_definitions(feast_mpi_obj PRIVATE MPI)
    add_library(pfeast $<TARGET_OBJECTS:feast_mpi_obj>)
endif ()

message(STATUS "FEAST Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")