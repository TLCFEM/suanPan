target_sources(${PROJECT_NAME} PRIVATE
        Constraint/Criterion/Criterion.cpp
        Constraint/Criterion/MaxDisplacement.cpp
        Constraint/Criterion/MaxHistory.cpp
        Constraint/Criterion/MaxResistance.cpp
        Constraint/Criterion/MinDisplacement.cpp
        Constraint/Criterion/MinResistance.cpp
        Constraint/Criterion/NodeBasedCriterion.cpp
        Constraint/Criterion/StrainEnergyEvolution.cpp
        Constraint/Criterion/ComplementaryEnergyEvolution.cpp
        Constraint/Criterion/EnergyEvolution.cpp
        )