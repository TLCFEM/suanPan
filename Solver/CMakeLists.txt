cmake_minimum_required(VERSION 3.13.0)

project(Solver)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../Libs)

include(Integrator/CMakeLists.txt)

add_library(${PROJECT_NAME} STATIC
        Arnoldi.cpp
        FEAST.cpp
        BFGS.cpp
        MPDC.cpp
        Newton.cpp
        Ramm.cpp
        Solver.cpp
        SolverParser.cpp
        ${Integrator}
        )

if (HAVE_VTK)
    target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})
    if (VTK_VERSION VERSION_LESS "8.90.0")
        include(${VTK_USE_FILE})
    else ()
        vtk_module_autoinit(TARGETS ${PROJECT_NAME} MODULES ${VTK_LIBRARIES})
    endif ()
endif ()
