target_sources(${PROJECT_NAME} PRIVATE
        Constraint.cpp
        ConstraintParser.cpp
        Embed.cpp
        FixedLength.cpp
        MPC.cpp
        NodeFacet.cpp
        NodeLine.cpp
        ParticleCollision.cpp
        ParticleCollision2D.cpp
        LJPotential2D.cpp
        LinearSpring2D.cpp
        ParticleCollision3D.cpp
        RestitutionWallPenalty.cpp
        RigidWallMultiplier.cpp
        RigidWallPenalty.cpp
)

add_subdirectory(BC)
add_subdirectory(Criterion)
