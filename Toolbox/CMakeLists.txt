target_sources(${PROJECT_NAME} PRIVATE
        argument.cpp
        command.cpp
        Expression.cpp
        ExpressionParser.cpp
        IntegrationPlan.cpp
        misc.cpp
        resampling.cpp
        sat.cpp
        sort_rcm.cpp
        tensor.cpp
        utility.cpp
)

if (SP_ENABLE_MPI)
    target_sources(${PROJECT_NAME} PRIVATE parpack.cpp)
else ()
    target_sources(${PROJECT_NAME} PRIVATE arpack.cpp)
endif ()
