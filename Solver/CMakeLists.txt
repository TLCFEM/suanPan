add_subdirectory(Integrator)

target_sources(${PROJECT_NAME} PRIVATE
        Arnoldi.cpp
        BFGS.cpp
        FEAST.cpp
        MPDC.cpp
        Newton.cpp
        Ramm.cpp
        Solver.cpp
        SolverParser.cpp
)
