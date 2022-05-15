target_sources(${PROJECT_NAME} PRIVATE
        BodyForce.cpp
        GroupBodyForce.cpp
        GroupNodalDisplacement.cpp
        GroupNodalForce.cpp
        Load.cpp
        NodalAcceleration.cpp
        NodalDisplacement.cpp
        NodalForce.cpp
        SupportMotion.cpp
        LineUDL.cpp
        )

add_subdirectory(Amplitude)