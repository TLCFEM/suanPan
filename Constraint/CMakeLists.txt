target_sources(${PROJECT_NAME} PRIVATE
        BC.cpp
        Constraint.cpp
        ConstraintParser.cpp
        Embed.cpp
        LinearSpring2D.cpp
        LJPotential2D.cpp
        MPC.cpp
        NodeFacet.cpp
        NodeLine.cpp
        ParticleCollision.cpp
        ParticleCollision2D.cpp
        ParticleCollision3D.cpp
        RestitutionWallPenalty.cpp
        RigidWallMultiplier.cpp
        RigidWallPenalty.cpp
)

add_subdirectory(Criterion)
