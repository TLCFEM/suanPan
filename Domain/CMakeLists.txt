target_sources(${PROJECT_NAME} PRIVATE
        Domain/Domain.cpp
        Domain/DomainBase.cpp
        Domain/ElementGroup.cpp
        Domain/ExternalModule.cpp
        Domain/Factory.hpp
        Domain/Group.cpp
        Domain/GroupGroup.cpp
        Domain/Node.cpp
        Domain/NodeGroup.cpp
        Domain/Storage.hpp
        Domain/Tag.cpp
        Domain/ConditionalModifier.cpp
        )

add_subdirectory(MetaMat)