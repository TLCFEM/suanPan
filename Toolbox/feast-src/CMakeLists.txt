cmake_minimum_required(VERSION 3.18)

project(feast Fortran)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

file(GLOB SPIKE_SRC CONFIGURE_DEPENDS
        spike-src/*.f
        spike-src/*.f90
)

set(FEAST_COMMON_SRC
        dzlsprim.f90
        feast_aux.f90
        sclsprim.f90
)

file(GLOB FEAST_SRC *.f90)
list(FILTER FEAST_SRC EXCLUDE REGEX "pdz")
list(FILTER FEAST_SRC EXCLUDE REGEX "(dzlsprim|feast_aux|sclsprim)\\.f90$")

file(GLOB FEAST_MPI_SRC *.f90)
list(FILTER FEAST_MPI_SRC EXCLUDE REGEX "(dzlsprim|feast_aux|sclsprim)\\.f90$")

if (SP_ENABLE_MKL)
    add_compile_definitions(MKL)
endif ()

if (NOT COMPILER_IDENTIFIER MATCHES "vs") # GNU GCC COMPILER
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -ffixed-line-length-none")
endif ()

add_library(feast_comm_obj OBJECT ${FEAST_COMMON_SRC} ${SPIKE_SRC})

if (SP_ENABLE_MPI)
    add_library(feast_mpi_obj OBJECT ${FEAST_MPI_SRC})
    target_compile_definitions(feast_mpi_obj PRIVATE MPI)
    target_link_libraries(feast_mpi_obj PUBLIC MPI::MPI_Fortran)

    add_library(pfeast $<TARGET_OBJECTS:feast_mpi_obj> $<TARGET_OBJECTS:feast_comm_obj>)
else ()
    add_library(feast_obj OBJECT ${FEAST_SRC})

    add_library(${PROJECT_NAME} $<TARGET_OBJECTS:feast_obj> $<TARGET_OBJECTS:feast_comm_obj>)
endif ()

message(STATUS "FEAST Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
