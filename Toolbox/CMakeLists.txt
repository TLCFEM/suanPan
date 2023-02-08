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
    target_compile_options(expression PRIVATE "/bigobj")
    if (CMAKE_BUILD_TYPE MATCHES Debug)
        target_compile_options(expression PRIVATE "/O1")
    endif ()
elseif (MINGW)
    target_compile_options(expression PRIVATE "-Wa,-mbig-obj")
    if (CMAKE_BUILD_TYPE MATCHES Debug)
        target_compile_options(expression PRIVATE "-O1")
    endif ()
endif ()

target_link_libraries(${PROJECT_NAME} expression)
