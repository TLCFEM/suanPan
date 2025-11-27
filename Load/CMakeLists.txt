target_sources(${PROJECT_NAME} PRIVATE
        BodyForce.cpp
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
