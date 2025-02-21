cmake_minimum_required(VERSION 3.18)

project(mpi-solver)

add_executable(solver.pardiso solver.pardiso.cpp)
add_executable(tester.pardiso tester.pardiso.cpp)

target_link_libraries(solver.pardiso ${MPI_LIBS})
target_link_libraries(tester.pardiso ${MPI_LIBS})

add_dependencies(tester.pardiso solver.pardiso)

add_dependencies(suanPan solver.pardiso)
