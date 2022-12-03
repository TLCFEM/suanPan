cmake_minimum_required(VERSION 3.13.0)

project(Converger)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../Libs)

add_library(${PROJECT_NAME} STATIC
        AbsDisp.cpp
        AbsError.cpp
        AbsIncreDisp.cpp
        AbsIncreAcc.cpp
        AbsIncreEnergy.cpp
        AbsResidual.cpp
        Converger.cpp
        ConvergerParser.cpp
        FixedNumber.cpp
        LogicConverger.cpp
        RelDisp.cpp
        RelError.cpp
        RelIncreDisp.cpp
        RelIncreAcc.cpp
        RelIncreEnergy.cpp
        RelResidual.cpp
        )
