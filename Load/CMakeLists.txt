target_sources(${PROJECT_NAME} PRIVATE
        BodyForce.cpp
        GroupBodyForce.cpp
        GroupNodalDisplacement.cpp
        GroupNodalForce.cpp
        LineUDL.cpp
        Load.cpp
        LoadParser.cpp
        NodalAcceleration.cpp
        NodalDisplacement.cpp
        NodalForce.cpp
        ReferenceForce.cpp
        SupportMotion.cpp
)

add_subdirectory(Amplitude)
