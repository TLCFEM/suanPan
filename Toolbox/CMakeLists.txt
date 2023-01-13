target_sources(${PROJECT_NAME} PRIVATE
        argument.cpp
        arpack.cpp
        command.cpp
        Converter.cpp
        IntegrationPlan.cpp
        sort_rcm.cpp
        sync_ostream.cpp
        tensor.cpp
        utility.cpp
        resampling.cpp
        )

add_library(expression STATIC Expression.cpp ExpressionParser.cpp)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/Libs)

if (MSVC)
    set_source_files_properties(Expression.cpp PROPERTY COMPILE_FLAGS /bigobj)
endif ()

target_link_libraries(${PROJECT_NAME} expression)
