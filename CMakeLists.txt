cmake_minimum_required(VERSION 3.0)

if(${CMAKE_VERSION} VERSION_GREATER "3.13.0")
    cmake_policy(SET CMP0076 OLD)
endif()

project(suanPan C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(CheckLanguage)
check_language(Fortran)
if(CMAKE_Fortran_COMPILER)
    set(FORTRAN_STATUS TRUE)
    enable_language(Fortran)
else()
    set(FORTRAN_STATUS FALSE)
endif()

set(ROOT .)

include(${ROOT}/.config.cmake)

link_directories(${PROJECT_BINARY_DIR}/Libs)

add_executable(${PROJECT_NAME}
        "suanPan.cpp"
        "Constraint/CMakeLists.txt"
        "Database/CMakeLists.txt"
        "Domain/CMakeLists.txt"
        "Load/CMakeLists.txt"
        "Recorder/CMakeLists.txt"
        "Step/CMakeLists.txt"
        "Toolbox/CMakeLists.txt"
        "UnitTest/CMakeLists.txt")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU") # GNU GCC COMPILER
    target_sources(${PROJECT_NAME} PRIVATE "Resource/suanPan_gcc.rc")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC") # MSVC COMPILER
    target_sources(${PROJECT_NAME} PRIVATE "Resource/suanPan.rc")
endif()

add_subdirectory(Constraint)
add_subdirectory(Converger)
add_subdirectory(Database)
add_subdirectory(Domain)
add_subdirectory(Element)
add_subdirectory(Load)
add_subdirectory(Material)
add_subdirectory(Recorder)
add_subdirectory(Section)
add_subdirectory(Solver)
add_subdirectory(Step)
add_subdirectory(Toolbox)
add_subdirectory(UnitTest)

if(FORTRAN_STATUS)
    add_subdirectory(Toolbox/amd-src)
    add_subdirectory(Toolbox/arpack-src)
    add_subdirectory(Toolbox/feast-src)
    add_subdirectory(Toolbox/lapack-ext)
    add_subdirectory(Toolbox/mumps-src)
endif()

target_link_libraries(${PROJECT_NAME}
    Converger
    Element
    Material
    Section
    Solver)

if(USE_SUPERLUMT)
    add_subdirectory(Toolbox/superlumt-src)
    target_link_libraries(${PROJECT_NAME} superlumt)
else()
    add_subdirectory(Toolbox/superlu-src)
    target_link_libraries(${PROJECT_NAME} superlu)
endif()

set_property(TARGET ${PROJECT_NAME} PROPERTY ENABLE_EXPORTS 1)
	
if((NOT FORTRAN_STATUS) AND (CMAKE_CXX_COMPILER_ID MATCHES "MSVC"))
    target_link_libraries(${PROJECT_NAME} libamd libarpack libspmm libmumps libfeast)
else()
    target_link_libraries(${PROJECT_NAME} amd arpack spmm mumps feast)
endif()

if((COMPILER_IDENTIFIER MATCHES "gcc-linux") OR (COMPILER_IDENTIFIER MATCHES "gcc-mac"))
    if(USE_MKL)
        set(DSUFIX ".a")
        set(DPREFIX "lib")
        if(BUILD_SHARED)
            set(DSUFIX "")
            set(DPREFIX "")
        endif()
        target_link_libraries(${PROJECT_NAME} ${DPREFIX}mkl_intel_lp64${DSUFIX} ${DPREFIX}mkl_intel_thread${DSUFIX} ${DPREFIX}mkl_core${DSUFIX} iomp5 m)
    else()
        target_link_libraries(${PROJECT_NAME} openblas)
    endif()
elseif(COMPILER_IDENTIFIER MATCHES "gcc-win")
    if(USE_MKL)
        target_link_libraries(${PROJECT_NAME} mkl_rt)
    else()
        target_link_libraries(${PROJECT_NAME} openblas)
    endif()
elseif(COMPILER_IDENTIFIER MATCHES "vs")
    if(USE_MKL)
        set(DSUFIX "")
        if(BUILD_SHARED)
            set(DSUFIX "_dll")
        endif()
        target_link_libraries(${PROJECT_NAME} mkl_intel_lp64${DSUFIX} mkl_intel_thread${DSUFIX} mkl_core${DSUFIX} libiomp5md)
    else()
        target_link_libraries(${PROJECT_NAME} libopenblas)
    endif()
endif()

if(BUILD_DLL_EXAMPLE)
    add_subdirectory(Developer/Element)
    add_subdirectory(Developer/Material)
    add_subdirectory(Developer/Modifier)
    add_subdirectory(Developer/ModuleBundle)
    add_subdirectory(Developer/Section)
    add_dependencies(ElementExample ${PROJECT_NAME})
    add_dependencies(MaterialExample ${PROJECT_NAME})
    add_dependencies(ModifierExample ${PROJECT_NAME})
    add_dependencies(ModuleBundle ${PROJECT_NAME})
    add_dependencies(SectionExample ${PROJECT_NAME})
endif()

# need further work
install(TARGETS ${PROJECT_NAME} DESTINATION ${ROOT}/bin)
install(FILES README.md CHANGELOG.md LICENSE DESTINATION ${ROOT}/bin)
install(FILES Enhancement/suanPan.sublime-completions Enhancement/suanPan.sublime-syntax DESTINATION ${ROOT}/bin)
if(COMPILER_IDENTIFIER MATCHES "gcc-linux")
    install(FILES Enhancement/suanPan.sh DESTINATION ${ROOT}/bin)
    install(FILES Libs/gcc-linux/libtbb.so.2 Libs/gcc-linux/libtbbmalloc.so.2 Libs/gcc-linux/libtbbmalloc_proxy.so.2 DESTINATION ${ROOT}/bin)
elseif(COMPILER_IDENTIFIER MATCHES "gcc-win")
install(FILES Enhancement/AddAssociation.bat DESTINATION ${ROOT}/bin)
    install(FILES Libs/gcc-win/tbb.dll Libs/gcc-win/tbbmalloc.dll Libs/gcc-win/tbbmalloc_proxy.dll DESTINATION ${ROOT}/bin)
elseif(COMPILER_IDENTIFIER MATCHES "vs")
    install(FILES Enhancement/AddAssociation.bat DESTINATION ${ROOT}/bin)
    install(FILES Libs/vs/tbb.dll Libs/vs/tbbmalloc.dll Libs/vs/tbbmalloc_proxy.dll DESTINATION ${ROOT}/bin)
endif()

if(HAVE_VTK)
    target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})
    if(VTK_VERSION VERSION_LESS "8.90.0")
        include(${VTK_USE_FILE})
    else()
        vtk_module_autoinit(TARGETS ${PROJECT_NAME} MODULES ${VTK_LIBRARIES})
    endif()
endif()

message("suanPan C_FLAGS: ${CMAKE_C_FLAGS}")
message("suanPan CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
if(FORTRAN_STATUS)
    message("suanPan Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
endif()

message("Link Dir:")
get_property(SDIR DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY LINK_DIRECTORIES)
foreach(SDIRA ${SDIR})
    message("${SDIRA}")
endforeach()

message("Include Dir:")
get_property(SDIR DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(SDIRA ${SDIR})
    message("${SDIRA}")
endforeach()