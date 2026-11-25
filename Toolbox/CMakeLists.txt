target_sources(${PROJECT_NAME} PRIVATE
        argument.cpp
        async_print.cpp
        command.cpp
        Expression.cpp
        ExpressionParser.cpp
        HDF5Wrapper.cpp
        IntegrationPlan.cpp
        misc.cpp
        resampling.cpp
        sort_rcm.cpp
        tensor.cpp
        utility.cpp
)

if (SP_ENABLE_MPI)
    target_sources(${PROJECT_NAME} PRIVATE eigen_parpack.cpp eigen_pfeast.cpp)
else ()
    target_sources(${PROJECT_NAME} PRIVATE eigen_arpack.cpp eigen_feast.cpp)
endif ()
