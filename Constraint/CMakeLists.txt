target_sources(${PROJECT_NAME} PRIVATE
        Constraint.cpp
        ConstraintParser.cpp
        Embed.cpp
        MPC.cpp
        ParticleCollision.cpp
        ParticleCollision2D.cpp
        ParticleCollision3D.cpp
        RigidWallMultiplier.cpp
        RigidWallPenalty.cpp
        FixedLength.cpp
        NodeLine.cpp
        NodeFacet.cpp
        )

add_subdirectory(BC)
add_subdirectory(Criterion)