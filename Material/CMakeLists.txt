add_subdirectory(Material1D)
add_subdirectory(Material2D)
add_subdirectory(Material3D)
add_subdirectory(MaterialOS)
add_subdirectory(Special)

target_sources(${PROJECT_NAME} PRIVATE
        ExternalMaterial.cpp
        Material.cpp
        MaterialParser.cpp
        MaterialTemplate.cpp
        MaterialTester.cpp
)
