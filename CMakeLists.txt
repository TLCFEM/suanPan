cmake_minimum_required(VERSION 3.18)

cmake_policy(SET CMP0076 NEW)

if (CMAKE_GNUtoMS_VCVARS OR MINGW)
    set(CMAKE_GNUtoMS ON)
endif ()

project(suanPan C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(CheckLanguage)
check_language(Fortran)
if (CMAKE_Fortran_COMPILER)
    set(FORTRAN_STATUS TRUE)
    enable_language(Fortran)
else ()
    set(FORTRAN_STATUS FALSE)
endif ()

set(ROOT ${CMAKE_CURRENT_SOURCE_DIR})

# run script to add revision tags to source file
if (WIN32)
    execute_process(COMMAND powershell ${ROOT}/Script/Rev.ps1 WORKING_DIRECTORY ${ROOT})
else ()
    execute_process(COMMAND bash -c ${ROOT}/Script/Rev.sh WORKING_DIRECTORY ${ROOT})
endif ()

# make sure changes to revision.h is not logged
execute_process(COMMAND git update-index --assume-unchanged Toolbox/revision.h WORKING_DIRECTORY ${ROOT})

include(${ROOT}/Option.cmake)

link_directories(${PROJECT_BINARY_DIR}/Libs)

add_executable(${PROJECT_NAME}
        suanPan.cpp
        Include/fmt/src/format.cc
        Include/whereami/whereami.c
        Include/catch/catch_amalgamated.cpp
        Constraint/CMakeLists.txt
        Database/CMakeLists.txt
        Domain/CMakeLists.txt
        Load/CMakeLists.txt
        Recorder/CMakeLists.txt
        Step/CMakeLists.txt
        Toolbox/CMakeLists.txt
        UnitTest/CMakeLists.txt
)

if (COMPILER_IDENTIFIER MATCHES "vs")
    target_link_options(${PROJECT_NAME} PRIVATE /NODEFAULTLIB:LIBCMT)
endif ()

if (MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE "/bigobj")
elseif (MINGW)
    target_compile_options(${PROJECT_NAME} PRIVATE "-Wa,-mbig-obj")
endif ()

if (COMPILER_IDENTIFIER MATCHES "vs") # MSVC COMPILER
    target_sources(${PROJECT_NAME} PRIVATE "Resource/suanPan.rc")
else ()
    target_sources(${PROJECT_NAME} PRIVATE "Resource/suanPan_gcc.rc")
endif ()

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

target_link_libraries(${PROJECT_NAME} Element Material Section Solver)

if (FORTRAN_STATUS)
    message(STATUS "Linking additional amd arpack mumps feast libraries.")
    add_subdirectory(Toolbox/amd-src)
    add_subdirectory(Toolbox/arpack-src)
    add_subdirectory(Toolbox/feast-src)
    add_subdirectory(Toolbox/mumps-src)
    add_subdirectory(Toolbox/fext)
    target_link_libraries(${PROJECT_NAME} fext)
elseif (COMPILER_IDENTIFIER MATCHES "vs")
    target_link_libraries(${PROJECT_NAME} libfext)
    message(STATUS "Linking precompiled fext (packed with amd arpack mumps feast) library.")
else ()
    message(FATAL_ERROR "Please install a valid FORTRAN compiler.")
endif ()

if (USE_SUPERLUMT)
    add_subdirectory(Toolbox/superlumt-src)
    target_link_libraries(${PROJECT_NAME} superlumt)
else ()
    add_subdirectory(Toolbox/superlu-src)
    target_link_libraries(${PROJECT_NAME} superlu)
endif ()

add_subdirectory(Toolbox/lis-src)
add_subdirectory(Toolbox/metis-src)
target_link_libraries(${PROJECT_NAME} lis metis)

set_property(TARGET ${PROJECT_NAME} PROPERTY ENABLE_EXPORTS 1)

if (USE_MPI)
    if (NOT COMPILER_IDENTIFIER MATCHES "(linux|mac)") # windows only support intel mpi
        set(MPI_LIBS impi impicxx)
    elseif (USE_INTEL_MPI)
        set(MPI_LIBS mpi mpifort)
    else ()
        set(MPI_LIBS mpi mpi_cxx mpi_mpifh)
    endif ()
    add_subdirectory(MPI)
    target_link_libraries(${PROJECT_NAME} ${MPI_LIBS})
endif ()

if (COMPILER_IDENTIFIER MATCHES "(linux|mac)")
    if (USE_MKL)
        if (USE_MPI)
            target_link_libraries(${PROJECT_NAME} MKL::MKL_SCALAPACK)
            target_link_libraries(solver.pardiso MKL::MKL_SCALAPACK)
            target_link_libraries(tester.pardiso MKL::MKL_SCALAPACK)
        else ()
            target_link_libraries(${PROJECT_NAME} MKL::MKL)
        endif ()
    elseif (USE_AOCL)
        if (NOT AOCL_BLIS_PATH MATCHES "blis")
            message(FATAL_ERROR "Please point AOCL_BLIS_PATH to the library file of BLIS.")
        endif ()
        if (NOT AOCL_FLAME_PATH MATCHES "flame")
            message(FATAL_ERROR "Please point AOCL_FLAME_PATH to the library file of FLAME.")
        endif ()
        if (NOT AOCL_UTILS_PATH MATCHES "aoclutils")
            message(FATAL_ERROR "Please point AOCL_UTILS_PATH to the library file of AOCL Utils.")
        endif ()
        target_link_libraries(${PROJECT_NAME} ${AOCL_FLAME_PATH} ${AOCL_BLIS_PATH} ${AOCL_UTILS_PATH})
    else ()
        target_link_libraries(${PROJECT_NAME} openblas)
        if (COMPILER_IDENTIFIER MATCHES "linux")
            set(CMAKE_EXE_LINKER_FLAGS "-Wl,-init,gotoblas_init")
        endif ()
    endif ()
elseif (COMPILER_IDENTIFIER MATCHES "win")
    if (USE_MKL)
        target_link_libraries(${PROJECT_NAME} mkl_rt)
    else ()
        target_link_libraries(${PROJECT_NAME} openblas)
    endif ()
elseif (COMPILER_IDENTIFIER MATCHES "vs")
    if (USE_MKL)
        set(DSUFIX "")
        if (LINK_DYNAMIC_MKL)
            set(DSUFIX "_dll")
        endif ()
        set(MKL_LIBS mkl_intel_lp64${DSUFIX} mkl_intel_thread${DSUFIX} mkl_core${DSUFIX} libiomp5md libifcore svml_disp)
        target_link_libraries(${PROJECT_NAME} ${MKL_LIBS})
        if (USE_MPI)
            if (USE_INTEL_MPI)
                list(APPEND MKL_LIBS mkl_blacs_intelmpi_lp64)
            else ()
                list(APPEND MKL_LIBS mkl_blacs_openmpi_lp64)
            endif ()
            target_link_libraries(solver.pardiso ${MKL_LIBS})
            target_link_libraries(tester.pardiso ${MKL_LIBS})
        endif ()
    else ()
        target_link_libraries(${PROJECT_NAME} libopenblas)
    endif ()
endif ()

if (BUILD_DLL_EXAMPLE)
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
endif ()

# need further work
install(TARGETS ${PROJECT_NAME} DESTINATION bin)
install(FILES
        CHANGELOG.md
        Enhancement/suanPan.sublime-completions
        Enhancement/suanPan.sublime-syntax
        LICENSE
        README.md
        Resource/suanPan-ua.svg
        Resource/suanPan.svg
        DESTINATION bin)
if (COMPILER_IDENTIFIER MATCHES "(linux|mac)")
    set(DECOR ".")
    set(SUFFIX "")
    if (COMPILER_IDENTIFIER MATCHES "linux")
        set(DECOR ".so")
        set(SUFFIX "")
    elseif (COMPILER_IDENTIFIER MATCHES "mac")
        set(DECOR "")
        set(SUFFIX ".dylib")
        file(GLOB OPENBLAS_FILES Libs/${SP_EXTERNAL_LIB_PATH}/libopenblas*)
        install(PROGRAMS ${OPENBLAS_FILES} DESTINATION lib)
    endif ()
    install(PROGRAMS Enhancement/suanPan.sh DESTINATION bin)
    if (BUILD_MULTITHREAD AND NOT USE_SYS_LIB)
        file(GLOB TBB_FILES Libs/${SP_EXTERNAL_LIB_PATH}/libtbb*)
        install(PROGRAMS ${TBB_FILES} DESTINATION lib)
    endif ()
    if (USE_MKL AND MKLROOT MATCHES "(oneapi|oneAPI)")
        if (LINK_DYNAMIC_MKL)
            file(GLOB MKL_FILES ${MKLROOT}/lib/intel64/libmkl_core${DECOR}*${SUFFIX})
            install(PROGRAMS ${MKL_FILES} DESTINATION lib)
            file(GLOB MKL_FILES ${MKLROOT}/lib/intel64/libmkl_intel_lp64${DECOR}*${SUFFIX})
            install(PROGRAMS ${MKL_FILES} DESTINATION lib)
            file(GLOB MKL_FILES ${MKLROOT}/lib/intel64/libmkl_def${DECOR}*${SUFFIX})
            install(PROGRAMS ${MKL_FILES} DESTINATION lib)
            file(GLOB MKL_FILES ${MKLROOT}/lib/intel64/libmkl_avx*${DECOR}*${SUFFIX})
            install(PROGRAMS ${MKL_FILES} DESTINATION lib)
            if (USE_INTEL_OPENMP)
                file(GLOB MKL_FILES ${MKLROOT}/lib/intel64/libmkl_intel_thread${DECOR}*${SUFFIX})
                install(PROGRAMS ${MKL_FILES} DESTINATION lib)
            else ()
                file(GLOB MKL_FILES ${MKLROOT}/lib/intel64/libmkl_gnu_thread${DECOR}*${SUFFIX})
                install(PROGRAMS ${MKL_FILES} DESTINATION lib)
            endif ()
        endif ()
        if (USE_INTEL_OPENMP)
            file(GLOB MKL_FILES ${IOMPPATH}/libiomp5${DECOR}*${SUFFIX})
            install(PROGRAMS ${MKL_FILES} DESTINATION lib)
        endif ()
    endif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
        set(PATH_LIST
                $ENV{ONEAPI_ROOT}/compiler/latest/lib #  since 2024
                $ENV{ONEAPI_ROOT}/compiler/latest/linux/compiler/lib/intel64_lin # prior to 2024
        )
        foreach (PATHA ${PATH_LIST})
            file(GLOB ONEAPI_FILES
                    ${PATHA}/libifcoremt${DECOR}*${SUFFIX}
                    ${PATHA}/libifport${DECOR}*${SUFFIX}
                    ${PATHA}/libimf${DECOR}*${SUFFIX}
                    ${PATHA}/libintlc${DECOR}*${SUFFIX}
                    ${PATHA}/libiomp5${DECOR}*${SUFFIX}
                    ${PATHA}/libsvml${DECOR}*${SUFFIX}
            )
            install(PROGRAMS ${ONEAPI_FILES} DESTINATION lib)
        endforeach ()
    endif ()
    if (USE_MAGMA)
        file(GLOB MAGMA_FILES ${MAGMAROOT}/lib/libmagma*.so)
        if (MAGMA_FILES)
            install(PROGRAMS ${MAGMA_FILES} DESTINATION lib)
        endif ()
    endif ()
elseif (COMPILER_IDENTIFIER MATCHES "win")
    install(FILES Enhancement/AddAssociation.bat DESTINATION bin)
    if (BUILD_MULTITHREAD)
        file(GLOB TBB_FILES Libs/${SP_EXTERNAL_LIB_PATH}/libtbb*.dll)
        install(FILES ${TBB_FILES} DESTINATION bin)
    endif ()
    if (USE_MKL AND MKLROOT MATCHES "(oneapi|oneAPI)")
        file(GLOB MKL_FILES
                ${MKLROOT}/redist/intel64/mkl_rt*
                ${MKLROOT}/redist/intel64/mkl_core*
                ${MKLROOT}/redist/intel64/mkl_def*
                ${MKLROOT}/redist/intel64/mkl_avx*
                ${MKLROOT}/redist/intel64/mkl_intel_thread.*
                ${MKLROOT}/bin/mkl_rt*
                ${MKLROOT}/bin/mkl_core*
                ${MKLROOT}/bin/mkl_def*
                ${MKLROOT}/bin/mkl_avx*
                ${MKLROOT}/bin/mkl_intel_thread.*
        )
        install(FILES ${MKL_FILES} DESTINATION bin)
    else ()
        file(GLOB DLL_FILES Libs/win/lib*.dll)
        install(FILES ${DLL_FILES} DESTINATION bin)
    endif ()
elseif (COMPILER_IDENTIFIER MATCHES "vs")
    install(FILES Enhancement/AddAssociation.bat DESTINATION bin)
    if (BUILD_MULTITHREAD OR USE_MKL)
        file(GLOB TBB_FILES Libs/${SP_EXTERNAL_LIB_PATH}/tbb*.dll)
        install(FILES ${TBB_FILES} DESTINATION bin)
    endif ()
    if (USE_MKL AND MKLROOT MATCHES "(oneapi|oneAPI)")
        if (LINK_DYNAMIC_MKL)
            file(GLOB MKL_FILES
                    ${MKLROOT}/redist/intel64/mkl_core*
                    ${MKLROOT}/redist/intel64/mkl_def*
                    ${MKLROOT}/redist/intel64/mkl_avx*
                    ${MKLROOT}/bin/mkl_core*
                    ${MKLROOT}/bin/mkl_def*
                    ${MKLROOT}/bin/mkl_avx*
            )
            install(FILES ${MKL_FILES} DESTINATION bin)
            file(GLOB MKL_FILES
                    ${MKLROOT}/redist/intel64/mkl_intel_thread.*
                    ${MKLROOT}/bin/mkl_intel_thread.*
            )
            install(FILES ${MKL_FILES} DESTINATION bin)
        endif ()
        if (USE_INTEL_OPENMP)
            find_file(libifcoremd libifcoremd.dll PATHS ${MKLROOT}/../../ ${MKLROOT}/../../compiler/latest/bin/ REQUIRED)
            install(FILES ${libifcoremd} DESTINATION bin)
            find_file(libiomp5md libiomp5md.dll PATHS ${MKLROOT}/../../ ${MKLROOT}/../../compiler/latest/bin/ REQUIRED)
            install(FILES ${libiomp5md} DESTINATION bin)
            find_file(libmmd libmmd.dll PATHS ${MKLROOT}/../../ ${MKLROOT}/../../compiler/latest/bin/ REQUIRED)
            install(FILES ${libmmd} DESTINATION bin)
            find_file(svml_dispmd svml_dispmd.dll PATHS ${MKLROOT}/../../ ${MKLROOT}/../../compiler/latest/bin/ REQUIRED)
            install(FILES ${svml_dispmd} DESTINATION bin)
        endif ()
    else ()
        file(GLOB DLL_FILES Libs/win/lib*.dll)
        install(FILES ${DLL_FILES} DESTINATION bin)
    endif ()
    file(GLOB VC_FILES Libs/vc/*.dll)
    install(FILES ${VC_FILES} DESTINATION bin)
endif ()

if (USE_MIMALLOC)
    add_dependencies(${PROJECT_NAME} mimalloc)
endif ()

if (HAVE_VTK)
    target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})
    if (VTK_VERSION VERSION_LESS "8.90.0")
        include(${VTK_USE_FILE})
    else ()
        vtk_module_autoinit(TARGETS ${PROJECT_NAME} MODULES ${VTK_LIBRARIES})
    endif ()
endif ()

message(STATUS "Flags and Dirs:")
message(STATUS "suanPan C_FLAGS: ${CMAKE_C_FLAGS}")
message(STATUS "suanPan CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
if (FORTRAN_STATUS)
    message(STATUS "suanPan Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
endif ()

message(STATUS "Link Dirs:")
get_property(SDIR DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY LINK_DIRECTORIES)
foreach (SDIRA ${SDIR})
    message(STATUS "${SDIRA}")
endforeach ()

message(STATUS "Include Dirs:")
get_property(SDIR DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach (SDIRA ${SDIR})
    message(STATUS "${SDIRA}")
endforeach ()

if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    file(READ "/etc/os-release" DISTRO_INFO)
    string(REGEX MATCH "fedora|ubuntu|debian" DIST ${DISTRO_INFO})

    if (DIST OR BUILD_PACKAGE)
        if ((DIST STREQUAL "ubuntu") OR (DIST STREQUAL "debian") OR (BUILD_PACKAGE MATCHES "DEB"))
            message(STATUS "Build DEB Package For Distribution: ${DIST}")
            set(CPACK_GENERATOR "DEB")
            # set(CPACK_DEBIAN_PACKAGE_DEPENDS "libstdc++6 libglvnd0 libgfortran5 libgomp1")
            set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
            set(CPACK_DEBIAN_COMPRESSION_TYPE "xz")
        elseif ((DIST STREQUAL "fedora") OR (BUILD_PACKAGE MATCHES "RPM"))
            message(STATUS "Build RPM Package For Distribution: ${DIST}")
            set(CPACK_GENERATOR "RPM")
            # set(CPACK_RPM_CHANGELOG_FILE ${ROOT}/CHANGELOG.md)
            set(CPACK_RPM_PACKAGE_AUTOREQPROV 1)
            set(CPACK_RPM_PACKAGE_LICENSE "GPL-3.0")
            set(CPACK_RPM_PACKAGE_REQUIRES "libstdc++ libglvnd libomp")
        endif ()

        set(CPACK_PACKAGE_CONTACT "Theodore Chang")
        set(CPACK_PACKAGE_CHECKSUM "SHA256")
        set(CPACK_PACKAGE_ICON ${ROOT}/Resource/suanPan-ua.svg)
        set(CPACK_PACKAGE_RELEASE 1)
        set(CPACK_PACKAGE_VENDOR "tlcfem")
        set(CPACK_PACKAGE_VERSION "3.7.0")
        set(CPACK_PACKAGE_DESCRIPTION "An Open Source, Parallel and Heterogeneous Finite Element Analysis Framework")
        set(CPACK_PACKAGE_HOMEPAGE_URL "https://github.com/TLCFEM/suanPan")

        set(CPACK_PACKAGE_FILE_NAME "${PROJECT_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_PACKAGE_RELEASE}.${CMAKE_SYSTEM_PROCESSOR}")

        include(CPack)
    endif ()
endif ()
