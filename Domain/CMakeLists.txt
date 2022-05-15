target_sources(${PROJECT_NAME} PRIVATE
        Domain.cpp
        DomainBase.cpp
        DomainHelper.cpp
        DomainState.cpp
        ElementGroup.cpp
        ExternalModule.cpp
        FactoryHelper.hpp
        Factory.hpp
        Group.cpp
        GroupGroup.cpp
        Node.cpp
        NodeGroup.cpp
        Storage.hpp
        Tag.cpp
        ConditionalModifier.cpp
        )

add_subdirectory(MetaMat)