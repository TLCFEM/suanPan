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
elseif (MINGW)
    if (CMAKE_BUILD_TYPE MATCHES Debug)
        set_source_files_properties(Expression.cpp PROPERTY COMPILE_FLAGS "-O1 -Wa,-mbig-obj")
    else ()
        set_source_files_properties(Expression.cpp PROPERTY COMPILE_FLAGS "-O2 -Wa,-mbig-obj")
    endif ()
endif ()

target_link_libraries(${PROJECT_NAME} expression)
