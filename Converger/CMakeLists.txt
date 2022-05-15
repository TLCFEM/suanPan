cmake_minimum_required(VERSION 3.13.0)

project(Converger)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../Libs)

add_library(${PROJECT_NAME} STATIC
        Converger.cpp
        AbsDisp.cpp
        AbsError.cpp
        AbsIncreDisp.cpp
        AbsIncreEnergy.cpp
        AbsResidual.cpp
        RelDisp.cpp
        RelError.cpp
        RelIncreDisp.cpp
        RelIncreEnergy.cpp
        RelResidual.cpp
        FixedNumber.cpp
        LogicConverger.cpp
        )