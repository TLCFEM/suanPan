add_subdirectory(Section1D)
add_subdirectory(Section2D)
add_subdirectory(Section3D)
add_subdirectory(SectionNM)
add_subdirectory(SectionOS)
add_subdirectory(SectionShell)

target_sources(${PROJECT_NAME} PRIVATE
        Fibre.cpp
        Section.cpp
        SectionParser.cpp
        SectionTester.cpp
)
