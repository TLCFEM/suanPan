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

if (USE_MKL)
    add_definitions(-DMKL)
endif ()

if (NOT COMPILER_IDENTIFIER MATCHES "vs") # GNU GCC COMPILER
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -ffixed-line-length-none")
endif ()

add_library(feast_obj OBJECT ${FEAST_SRC} ${SPIKE_SRC})

add_library(${PROJECT_NAME} ${LIBRARY_TYPE} $<TARGET_OBJECTS:feast_obj>)

message(STATUS "FEAST Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")