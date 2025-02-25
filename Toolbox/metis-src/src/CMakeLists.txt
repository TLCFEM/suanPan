include_directories(.)

file(GLOB METIS_SRC *.c)

add_library(${PROJECT_NAME} ${GK_SRC} ${METIS_SRC})

message(STATUS "METIS C_FLAGS: ${CMAKE_C_FLAGS}")

if (UNIX)
    target_link_libraries(${PROJECT_NAME} m)
endif ()
