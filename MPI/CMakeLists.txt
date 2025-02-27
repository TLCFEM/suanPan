cmake_minimum_required(VERSION 3.18)

project(mpi-solver)

if (CMAKE_BUILD_TYPE MATCHES "Deb")
    add_compile_definitions(MPL_DEBUG)
endif()

add_executable(solver.pardiso solver.pardiso.cpp)
add_executable(tester.pardiso tester.pardiso.cpp)

target_link_libraries(solver.pardiso MKL::MKL_SCALAPACK MPI::MPI_C MPI::MPI_CXX MPI::MPI_Fortran)
target_link_libraries(tester.pardiso MKL::MKL_SCALAPACK MPI::MPI_C MPI::MPI_CXX MPI::MPI_Fortran)

add_dependencies(tester.pardiso solver.pardiso)

add_dependencies(suanPan solver.pardiso)

install(TARGETS solver.pardiso DESTINATION bin)
