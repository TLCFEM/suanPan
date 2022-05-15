cmake_minimum_required(VERSION 3.13.0)

project(Element)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../Libs)

include(Beam/CMakeLists.txt)
include(Cube/CMakeLists.txt)
include(Membrane/CMakeLists.txt)
include(Modifier/CMakeLists.txt)
include(Plate/CMakeLists.txt)
include(Shell/CMakeLists.txt)
include(Special/CMakeLists.txt)
include(Truss/CMakeLists.txt)
include(Utility/CMakeLists.txt)
include(Visualisation/CMakeLists.txt)
include(Patch/CMakeLists.txt)

add_library(${PROJECT_NAME} STATIC
        Element.cpp
        ElementParser.cpp
        ElementTemplate.cpp
        MaterialElement.cpp
        SectionElement.cpp
        ${Beam}
        ${Cube}
        ${Membrane}
        ${Modifier}
        ${Plate}
        ${Shell}
        ${Special}
        ${Truss}
        ${Utility}
        ${Visualisation}
        ${Patch}
        )

if (HAVE_VTK)
    target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})
    if (VTK_VERSION VERSION_LESS "8.90.0")
        include(${VTK_USE_FILE})
    else ()
        vtk_module_autoinit(TARGETS ${PROJECT_NAME} MODULES ${VTK_LIBRARIES})
    endif ()
endif ()