add_subdirectory(Beam)
add_subdirectory(Cube)
add_subdirectory(Membrane)
add_subdirectory(Modifier)
add_subdirectory(Patch)
add_subdirectory(Plate)
add_subdirectory(Shell)
add_subdirectory(Special)
add_subdirectory(Truss)
add_subdirectory(Utility)
add_subdirectory(Visualisation)

target_sources(${PROJECT_NAME} PRIVATE
        Element.cpp
        ElementParser.cpp
        ElementTemplate.cpp
        MaterialElement.cpp
        SectionElement.cpp
)
