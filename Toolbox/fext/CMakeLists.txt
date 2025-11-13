cmake_minimum_required(VERSION 3.18)

project(fext Fortran)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

add_library(${PROJECT_NAME}
        $<TARGET_OBJECTS:arpack_obj>
        $<TARGET_OBJECTS:feast_obj>
)

if (BUILD_SHARED_LIBS)
    if (WIN32)
        install(TARGETS ${PROJECT_NAME} DESTINATION bin COMPONENT Main)
    else ()
        install(TARGETS ${PROJECT_NAME} DESTINATION lib COMPONENT Main)
    endif ()
endif ()
