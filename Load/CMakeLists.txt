target_sources(${PROJECT_NAME} PRIVATE
        "Load/BodyForce.cpp"
        "Load/GroupBodyForce.cpp"
        "Load/GroupNodalDisplacement.cpp"
        "Load/GroupNodalForce.cpp"
        "Load/Load.cpp"
        "Load/NodalAcceleration.cpp"
        "Load/NodalDisplacement.cpp"
        "Load/NodalForce.cpp"
        "Load/SupportMotion.cpp"
        )

add_subdirectory(Amplitude)