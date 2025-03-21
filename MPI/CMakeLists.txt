cmake_minimum_required(VERSION 3.18)

project(mpi-solver)

if (CMAKE_BUILD_TYPE MATCHES "Deb")
    add_compile_definitions(MPL_DEBUG)
endif()

add_executable(distributed_obj distributed_obj.cpp)
target_link_libraries(distributed_obj MKL::MKL_SCALAPACK MPI::MPI_C MPI::MPI_CXX MPI::MPI_Fortran)
