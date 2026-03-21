cmake_minimum_required(VERSION 3.18)

cmake_policy(SET CMP0076 NEW)

project(suanPan C CXX)

if (CMAKE_GNUtoMS_VCVARS AND MINGW)
    set(CMAKE_GNUtoMS ON)
endif ()

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(CheckLanguage)
check_language(Fortran)
if (CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
endif ()

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
    # make sure changes to revision.h is not logged
    execute_process(COMMAND git update-index --assume-unchanged Toolbox/revision.h WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
    find_program(GIT_EXE NAMES git)
    if (GIT_EXE)
        execute_process(
                COMMAND git rev-parse --short=8 HEAD
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                OUTPUT_VARIABLE GIT_COMMIT_HASH
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
        )
        if (GIT_COMMIT_HASH)
            configure_file(
                    ${CMAKE_CURRENT_SOURCE_DIR}/Toolbox/revision.h.in
                    ${CMAKE_CURRENT_SOURCE_DIR}/Toolbox/revision.h
                    @ONLY
            )
            message(STATUS "Revision set to ${GIT_COMMIT_HASH}.")
        endif ()
    else ()
        message(STATUS "Git executable not found.")
    endif ()
endif ()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(Include)
include_directories(Include/fmt/include)

message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

set(BUILD_PACKAGE "" CACHE STRING "DEB OR RPM")

if (CMAKE_SYSTEM_PROCESSOR MATCHES "x86" OR CMAKE_SYSTEM_PROCESSOR MATCHES "AMD64")
    set(SP_AMD64 ON)
endif ()

option(BUILD_SHARED_LIBS "Build shared instead of static libraries." OFF)
option(SP_BUILD_DLL_EXAMPLE "Build dynamic linked library examples." OFF)
option(SP_BUILD_PARALLEL "Build with parallel support via TBB." OFF)
option(SP_ENABLE_AVX "Build with AVX support, has no effect on ARM." OFF)
option(SP_ENABLE_AVX2 "Build with AVX2 support, has no effect on ARM." ON)
option(SP_ENABLE_AVX512 "Build with AVX512 support, has no effect on ARM." OFF)
option(SP_ENABLE_NATIVE_ARCH "Build with native arch, when this is set, all AVX related settings will be ignored." OFF)
option(SP_ENABLE_HDF5 "Build with HDF5 format support." ON)
option(SP_ENABLE_LTO "Enable link-time optimisation." OFF)
if (SP_AMD64)
    option(SP_ENABLE_MKL "Build with Intel MKL instead of the bundled OpenBLAS." OFF)
else ()
    unset(SP_ENABLE_MKL CACHE)
endif ()
option(SP_ENABLE_VTK "Build with visualisation support via VTK, note external VTK libraries need to be present." OFF)
if (LINUX AND CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    option(SP_ENABLE_AOCL "Build with AMD Optimizing CPU Libraries instead of the bundled OpenBLAS." OFF)
endif ()
option(SP_ENABLE_MIMALLOC "Build with mimalloc instead of the default memory allocator." OFF)
option(SP_ENABLE_SYSLIB "Build with libraries (HDF5 and TBB) installed on the system instead of the bundled ones." OFF)
if (SP_ENABLE_MKL)
    option(SP_ENABLE_64BIT_INDEXING "Use 64-bit integer for indexing." OFF)
    option(SP_ENABLE_IOMP "Use Intel's OpenMP implementation." OFF)
    option(SP_ENABLE_SHARED_MKL "Use dynamic Intel MKL libraries instead of static ones." OFF)
    option(SP_ENABLE_MPI "Build with cluster parallelism support via MPI, note a working MPI implementation and ScaLAPACK must be present." OFF)
else ()
    # set(SP_ENABLE_MPI OFF CACHE BOOL "" FORCE)
    set(SP_OPENBLAS_PATH "" CACHE FILEPATH "Manually assign the path to the OpenBLAS library, if not set, the bundled version will be used.")
endif ()
if (SP_ENABLE_AOCL)
    set(SP_AOCL_PATH "${CMAKE_SOURCE_DIR}/Libs/linux/libflame.a;${CMAKE_SOURCE_DIR}/Libs/linux/libblis-mt.a;${CMAKE_SOURCE_DIR}/Libs/linux/libaoclutils.a" CACHE FILEPATH "AOCL all-in-one library path.")
    add_compile_definitions(SUANPAN_AOCL)
endif ()

if (SUANPAN_DATE)
    add_compile_definitions(SUANPAN_DATE=${SUANPAN_DATE})
endif ()

set(COMPILER_IDENTIFIER "unknown")
set(SP_EXTERNAL_LIB_PATH "unknown")
if (CMAKE_SYSTEM_NAME MATCHES "Windows") # WINDOWS PLATFORM
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU") # GNU GCC COMPILER
        set(COMPILER_IDENTIFIER "gcc-win")
        set(SP_EXTERNAL_LIB_PATH "win")
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC" OR CMAKE_CXX_COMPILER_ID MATCHES "Intel") # MSVC COMPILER
        set(COMPILER_IDENTIFIER "vs")
        set(SP_EXTERNAL_LIB_PATH "vs")
        add_compile_definitions(_SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING)
        option(SP_ENABLE_CUDA "Enable GPU based global solvers via CUDA." OFF)
    endif ()
elseif (CMAKE_SYSTEM_NAME MATCHES "Linux") # LINUX PLATFORM
    set(SP_EXTERNAL_LIB_PATH "linux")
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU") # GNU GCC COMPILER
        set(COMPILER_IDENTIFIER "gcc-linux")
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM") # Intel COMPILER icpx
        set(COMPILER_IDENTIFIER "clang-linux")
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang") # Clang COMPILER
        set(COMPILER_IDENTIFIER "clang-linux")
    endif ()
    option(SP_ENABLE_CUDA "Enable GPU based global solvers via CUDA." OFF)
elseif (CMAKE_SYSTEM_NAME MATCHES "Darwin") # MAC PLATFORM
    set(SP_EXTERNAL_LIB_PATH "mac")
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU") # GNU GCC COMPILER
        set(COMPILER_IDENTIFIER "gcc-mac")
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(COMPILER_IDENTIFIER "clang-mac")
        include_directories(
                "/usr/local/include"
                "/usr/local/opt/llvm/include"
                "/usr/local/opt/libomp/include"
                "/opt/homebrew/opt/libomp/include"
        )
        link_directories(
                "/usr/local/lib"
                "/usr/local/opt/llvm/lib"
                "/usr/local/opt/libomp/lib/"
                "/opt/homebrew/opt/libomp/lib"
        )
        message(NOTICE "On macOS, make sure `libomp` is installed.")
    endif ()
endif ()

if (COMPILER_IDENTIFIER MATCHES "unknown")
    message(FATAL_ERROR "Cannot identify the compiler available, please use GCC or MSVC or Intel.")
endif ()

if (NOT SP_AMD64)
    set(SP_EXTERNAL_LIB_PATH "${SP_EXTERNAL_LIB_PATH}-arm")
    if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Libs/${SP_EXTERNAL_LIB_PATH}")
        message(NOTICE "Dependencies are not bundled and not downloaded, please check `Script/mac-dependency.sh` or `Script/linux-dependency.sh` for more details.")
    endif ()
endif ()

link_directories(Libs/${SP_EXTERNAL_LIB_PATH})

if (SP_ENABLE_MKL)
    if (SP_ENABLE_CUDA)
        option(SP_ENABLE_MAGMA "Enable GPU based global solvers via MAGMA." OFF)
    endif ()
    set(MKL_PATH "" CACHE PATH "MKL library path which contains /include and /lib folders.")
    if (CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
        set(SP_ENABLE_IOMP ON CACHE BOOL "" FORCE)
    endif ()
    if (COMPILER_IDENTIFIER MATCHES "win")
        set(MKL_LINK sdl)
    elseif (NOT SP_ENABLE_SHARED_MKL)
        set(MKL_LINK static)
    endif ()
    if (NOT SP_ENABLE_IOMP AND LINUX)
        set(MKL_THREADING gnu_thread)
    endif ()
    if (SP_ENABLE_64BIT_INDEXING)
        message(NOTICE "64-bit indexing is enabled.")
        add_compile_definitions(SUANSPAN_64BIT_INT)
        add_compile_definitions(EZP_INT64)
        set(MKL_INTERFACE ilp64)
    else ()
        set(MKL_INTERFACE lp64)
    endif ()
    if (SP_ENABLE_MPI)
        add_compile_definitions(SUANPAN_DISTRIBUTED MPICH_SKIP_MPICXX)
        find_package(MPI REQUIRED COMPONENTS C Fortran)
        get_filename_component(MPI_ROOT_DIR ${MPIEXEC_EXECUTABLE} DIRECTORY)
        file(GLOB MPI_VERSION_FILE
                "${MPI_ROOT_DIR}/mpichversion"
                "${MPI_ROOT_DIR}/ompi_info"
                "${MPI_ROOT_DIR}/impi_info"
        )
        if (NOT MKL_MPI)
            if (MPI_VERSION_FILE MATCHES "ompi_info")
                set(MKL_MPI openmpi)
            elseif (MPI_VERSION_FILE MATCHES "mpichversion")
                set(MKL_MPI mpich)
            endif ()
        endif ()
        set(EZP_ENABLE_OPENMP ON)
        add_subdirectory(Include/ezp)
        include_directories(Include/ezp/mpl)
    endif ()
    find_package(MKL REQUIRED PATHS ${MKL_PATH} "/opt/intel/oneapi/mkl/latest" "C:/Program Files (x86)/Intel/oneAPI/mkl/latest")
    add_compile_definitions(SUANPAN_MKL)
    add_compile_definitions(EZP_MKL)
    # add_compile_definitions(ARMA_USE_MKL_ALLOC)
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    message(WARNING "Since Intel compilers are used, why not enabling MKL?")
elseif (SP_ENABLE_64BIT_INDEXING)
    message(NOTICE "64-bit indexing is enabled, please make sure the linear algebra driver is compiled using 64-bit integer as well.")
endif ()

if (SP_ENABLE_CUDA)
    find_package(CUDAToolkit REQUIRED)
    link_libraries(CUDA::cudart_static CUDA::cublas_static CUDA::cusolver_static CUDA::cusparse_static)
    add_compile_definitions(SUANPAN_CUDA)
    if (SP_ENABLE_MAGMA)
        set(MAGMA_PATH "" CACHE PATH "Magma library path which contains /include and /lib folders.")
        find_file(MAGMA_HEADER NAMES magma.h PATHS ${MAGMA_PATH}/include NO_DEFAULT_PATH)
        if (MAGMA_HEADER MATCHES "MAGMA_HEADER-NOTFOUND")
            message(FATAL_ERROR "The MAGMA header <magma.h> is not found: ${MAGMA_PATH}.")
        endif ()
        include_directories(${MAGMA_PATH}/include)
        link_directories(${MAGMA_PATH}/lib)
        link_libraries(magma magma_sparse)
        add_compile_definitions(SUANPAN_MAGMA)
    endif ()
endif ()

if (SP_ENABLE_VTK)
    set(VTK_PATH "" CACHE PATH "VTK library path which contains /include folder.")
    if (SP_ENABLE_SYSLIB OR NOT VTK_PATH)
        find_package(VTK REQUIRED)
    else ()
        find_package(VTK REQUIRED PATHS ${VTK_PATH} NO_DEFAULT_PATH)
    endif ()
    add_compile_definitions(SUANPAN_VTK)
endif ()

if (SP_ENABLE_HDF5)
    add_compile_definitions(SUANPAN_HDF5)
    if (SP_ENABLE_VTK)
        string(REGEX REPLACE "/lib6?4?/cmake/vtk" "/include/vtk" VTK_INCLUDE ${VTK_DIR}) # on linux
        string(REGEX REPLACE "\\\\lib6?4?\\\\cmake\\\\vtk" "\\\\include\\\\vtk" VTK_INCLUDE ${VTK_INCLUDE}) # on windows
        include_directories(${VTK_INCLUDE}/vtkhdf5 ${VTK_INCLUDE}/vtkhdf5/src ${VTK_INCLUDE}/vtkhdf5/hl/src)
        find_file(HDF5_HEADER NAMES hdf5.h PATHS ${VTK_INCLUDE}/vtkhdf5/src NO_DEFAULT_PATH)
        if (HDF5_HEADER MATCHES "HDF5_HEADER-NOTFOUND")
            message(FATAL_ERROR "The HDF5 header <hdf5.h> is not found.")
        else ()
            message(STATUS "Found HDF5 header <hdf5.h>: ${HDF5_HEADER}")
        endif ()
    elseif (SP_ENABLE_SYSLIB)
        find_package(HDF5 REQUIRED COMPONENTS C HL)
        link_libraries(hdf5::hdf5 hdf5::hdf5_hl)
    else ()
        include_directories(Include/hdf5 Include/hdf5-${SP_EXTERNAL_LIB_PATH})
        if (COMPILER_IDENTIFIER MATCHES "vs")
            link_libraries(libhdf5_hl libhdf5 shlwapi)
        else ()
            link_libraries(hdf5_hl hdf5)
        endif ()
    endif ()
endif ()

if (SP_BUILD_PARALLEL)
    message(STATUS "USING TBB LIBRARY")
    add_compile_definitions(SUANPAN_MT)
    if (SP_ENABLE_SYSLIB)
        find_package(TBB REQUIRED)
        link_libraries(TBB::tbb TBB::tbbmalloc TBB::tbbmalloc_proxy)
    else ()
        include_directories(Include/tbb)
        if (COMPILER_IDENTIFIER MATCHES "gcc-win")
            link_libraries(tbb12)
        else ()
            link_libraries(tbb)
        endif ()
        option(SP_ENABLE_TBB_ALLOC "Build with TBB's memory allocator." OFF)
        if (SP_ENABLE_TBB_ALLOC)
            add_compile_definitions(ARMA_USE_TBB_ALLOC) # for armadillo to use tbb allocator
            include_directories(Include/tbb/oneapi) # because armadillo assumes oneapi be in the include path
            link_libraries(tbbmalloc tbbmalloc_proxy)
        endif ()
    endif ()
endif ()

if (CMAKE_BUILD_TYPE MATCHES "Rel")
    add_compile_definitions(ARMA_WARN_LEVEL=1 ARMA_DONT_CHECK_CONFORMANCE)
endif ()

if (BUILD_SHARED_LIBS)
    message(STATUS "BUILD SHARED LIBRARY")
else ()
    message(STATUS "BUILD STATIC LIBRARY")
endif ()

if (SP_AMD64)
    if (SP_ENABLE_AVX512)
        add_compile_definitions(SUANPAN_AVX512)
    elseif (SP_ENABLE_AVX2)
        add_compile_definitions(SUANPAN_AVX2)
    elseif (SP_ENABLE_AVX)
        add_compile_definitions(SUANPAN_AVX)
    endif ()
endif ()

if (SP_ENABLE_SUPERLUMT)
    message(WARNING "The multithreaded version of SuperLU tends to be slow, it is not recommended to use it.")
    add_compile_definitions(SUANPAN_SUPERLUMT)
endif ()

if (COMPILER_IDENTIFIER MATCHES "vs")
    unset(SP_ENABLE_CODECOV CACHE)

    link_directories(Libs/win)

    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MP")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /EHsc /wd4068 /wd4996")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /nowarn /MP /Qparallel /fpp /names:lowercase /assume:underscore /libs:dll /threads")
    if (SP_ENABLE_64BIT_INDEXING)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /integer-size:64")
    endif ()

    if (SP_AMD64)
        if (SP_ENABLE_NATIVE_ARCH)
            # do nothing
        elseif (SP_ENABLE_AVX512)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /arch:AVX512")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX512")
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /arch:AVX2")
        elseif (SP_ENABLE_AVX2)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /arch:AVX2")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2")
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /arch:AVX2")
        elseif (SP_ENABLE_AVX)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /arch:AVX")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX")
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /arch:AVX")
        endif ()
    endif ()
else ()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fexceptions")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fexceptions")
    if (CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
        link_libraries(stdc++)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffp-model=precise")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffp-model=precise")
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Werror -Wno-unused-parameter -Wno-parentheses -Wno-maybe-uninitialized")
    endif ()

    if (SP_AMD64)
        if (SP_ENABLE_NATIVE_ARCH)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
        elseif (SP_ENABLE_AVX512)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512f")
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx512f")
        elseif (SP_ENABLE_AVX2)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx2")
        elseif (SP_ENABLE_AVX)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx")
        endif ()
    endif ()

    option(SP_ENABLE_CODECOV "Enable code coverage report." OFF)

    if (SP_ENABLE_CODECOV AND COMPILER_IDENTIFIER MATCHES "gcc") # only report coverage with gcc
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage -fprofile-update=atomic")
        link_libraries(gcov)
    endif ()

    if (CMAKE_BUILD_TYPE MATCHES "Debug")
        option(SP_ENABLE_ASAN "Enable address sanitizer." OFF)
        if (SP_ENABLE_ASAN)
            message(STATUS "Using the address sanitizer with flags: -fsanitize=address,leak,undefined,bounds-strict,pointer-overflow,shift -fsanitize-address-use-after-scope -fno-omit-frame-pointer -fno-optimize-sibling-calls")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address,leak,undefined,bounds-strict,pointer-overflow,shift -fsanitize-address-use-after-scope -fno-omit-frame-pointer -fno-optimize-sibling-calls")
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address,leak,undefined,bounds-strict,pointer-overflow,shift -fsanitize-address-use-after-scope -fno-omit-frame-pointer -fno-optimize-sibling-calls")
        endif ()
    endif ()

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w -fallow-argument-mismatch")
    if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp")
        if (SP_ENABLE_64BIT_INDEXING)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -integer-size 64")
        endif ()

        if (SP_AMD64)
            if (SP_ENABLE_NATIVE_ARCH)
                # do nothing
            elseif (SP_ENABLE_AVX512)
                set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -arch AVX2")
            elseif (SP_ENABLE_AVX2)
                set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -arch AVX2")
            elseif (SP_ENABLE_AVX)
                set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -arch AVX")
            endif ()
        endif ()
    else ()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
        if (SP_ENABLE_64BIT_INDEXING)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-integer-8")
        endif ()

        if (SP_AMD64)
            if (SP_ENABLE_NATIVE_ARCH)
                set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=native")
            elseif (SP_ENABLE_AVX512)
                set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mavx512f")
            elseif (SP_ENABLE_AVX2)
                set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mavx2")
            elseif (SP_ENABLE_AVX)
                set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mavx")
            endif ()
        endif ()
    endif ()
endif ()

if (SP_ENABLE_MKL)
    if (SP_ENABLE_MPI)
        link_libraries(MKL::MKL_SCALAPACK lis mumps)
    else ()
        link_libraries(MKL::MKL)
    endif ()
elseif (SP_ENABLE_AOCL)
    link_libraries(${SP_AOCL_PATH})
else ()
    if (SP_OPENBLAS_PATH)
        link_libraries(${SP_OPENBLAS_PATH})
    elseif (SP_ENABLE_SYSLIB)
        find_library(OPENBLAS_LIB NAMES openblaso openblas)
        if (NOT OPENBLAS_LIB)
            message(FATAL_ERROR "Cannot find OpenBLAS library on the system. Alternatively, explicit assign the path via SP_OPENBLAS_PATH.")
        endif ()
        link_libraries(${OPENBLAS_LIB})
    elseif (COMPILER_IDENTIFIER MATCHES "vs")
        link_libraries(libopenblas)
    else ()
        link_libraries(openblas)
    endif ()
    if (CMAKE_SYSTEM_NAME MATCHES "Linux")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-init,gotoblas_init")
    endif ()
endif ()

if (MINGW)
    foreach (STATIC_DEP gfortran quadmath gomp dl)
        find_library(${STATIC_DEP}_LIB NAMES lib${STATIC_DEP}.a PATHS
                ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}
                ${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES}
                ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}
        )
        if (${STATIC_DEP}_LIB)
            message(STATUS "Found static library ${${STATIC_DEP}_LIB}")
            link_libraries(${${STATIC_DEP}_LIB})
        endif ()
    endforeach ()
endif ()

add_executable(${PROJECT_NAME}
        suanPan.cpp
        Include/catch/catch_amalgamated.cpp
        Include/fmt/src/format.cc
        Include/whereami/whereami.c
)

set_property(TARGET ${PROJECT_NAME} PROPERTY ENABLE_EXPORTS 1)

if ((NOT MINGW AND CMAKE_BUILD_TYPE STREQUAL "Release") OR SP_ENABLE_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT SP_LTO_RESULT OUTPUT SP_LTO_OUTPUT)
    if (SP_LTO_RESULT)
        set_property(TARGET ${PROJECT_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif ()
endif ()

if (MSVC)
    target_sources(${PROJECT_NAME} PRIVATE "Resource/suanPan.rc")
    target_compile_options(${PROJECT_NAME} PRIVATE /bigobj /utf-8)
    target_link_options(${PROJECT_NAME} PRIVATE "/NODEFAULTLIB:LIBCMT")
else ()
    target_sources(${PROJECT_NAME} PRIVATE "Resource/suanPan_gcc.rc")
    if (MINGW)
        target_compile_options(${PROJECT_NAME} PRIVATE "-Wa,-mbig-obj")
    endif ()
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

if (CMAKE_Fortran_COMPILER)
    find_package(OpenMP REQUIRED COMPONENTS Fortran)
    link_libraries(OpenMP::OpenMP_Fortran)
    message(STATUS "Linking additional arpack feast libraries.")
    add_subdirectory(Toolbox/arpack-src)
    add_subdirectory(Toolbox/feast-src)
    add_subdirectory(Toolbox/fext)
    target_link_libraries(${PROJECT_NAME} fext)
    if (SP_ENABLE_MPI)
        add_subdirectory(Toolbox/parpack-src)
        target_link_libraries(${PROJECT_NAME} parpack)
    endif ()
elseif (COMPILER_IDENTIFIER MATCHES "vs")
    target_link_libraries(${PROJECT_NAME} libfext)
    message(STATUS "Linking precompiled fext (packed with arpack feast) library.")
else ()
    message(FATAL_ERROR "Please install a valid FORTRAN compiler.")
endif ()

if (SP_ENABLE_SUPERLUMT)
    add_subdirectory(Toolbox/superlumt-src)
    target_link_libraries(${PROJECT_NAME} superlumt)
else ()
    add_subdirectory(Toolbox/superlu-src)
    target_link_libraries(${PROJECT_NAME} superlu)
endif ()

if (SP_ENABLE_MIMALLOC)
    message(STATUS "USING MIMALLOC LIBRARY")
    include(FetchContent)
    FetchContent_Declare(mimalloc
            GIT_REPOSITORY https://github.com/microsoft/mimalloc
            GIT_TAG v2.2.3)
    FetchContent_MakeAvailable(mimalloc)
    set(MI_BUILD_OBJECT OFF)
    set(MI_BUILD_STATIC OFF)
    set(MI_BUILD_TESTS OFF)
    set(MI_INSTALL_TOPLEVEL ON)
    add_dependencies(${PROJECT_NAME} mimalloc)
endif ()

if (SP_ENABLE_VTK)
    target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})
    if (VTK_VERSION VERSION_LESS "8.90.0")
        include(${VTK_USE_FILE})
    else ()
        vtk_module_autoinit(TARGETS ${PROJECT_NAME} MODULES ${VTK_LIBRARIES})
    endif ()
endif ()

if (SP_BUILD_DLL_EXAMPLE)
    add_subdirectory(Developer/Element)
    add_subdirectory(Developer/Material)
    add_subdirectory(Developer/Modifier)
    add_subdirectory(Developer/ModuleBundle)
    add_subdirectory(Developer/Section)
endif ()

if (MINGW AND NOT BUILD_SHARED_LIBS AND NOT SP_BUILD_PARALLEL)
    target_link_options(${PROJECT_NAME} PRIVATE -static)
endif ()

# need further work
install(TARGETS ${PROJECT_NAME} DESTINATION bin COMPONENT Main)
if (COMPILER_IDENTIFIER MATCHES "(linux|mac)")
    install(FILES
            Resource/suanPan-ua.svg
            DESTINATION share/icons/hicolor/scalable/apps
            RENAME suanPan.svg
            COMPONENT Support
    )
    install(FILES
            Enhancement/suanPan.sublime-completions
            Enhancement/suanPan.sublime-syntax
            DESTINATION share/${PROJECT_NAME}
            COMPONENT Support
    )
    install(FILES
            CHANGELOG.md
            LICENSE
            README.md
            DESTINATION share/doc/${PROJECT_NAME}
            COMPONENT Support
    )
else ()
    install(FILES
            CHANGELOG.md
            Enhancement/AddAssociation.bat
            Enhancement/suanPan.sublime-completions
            Enhancement/suanPan.sublime-syntax
            LICENSE
            README.md
            Resource/suanPan-ua.svg
            DESTINATION bin
            COMPONENT Support
    )
endif ()
if (COMPILER_IDENTIFIER MATCHES "(linux|mac)")
    install(PROGRAMS Enhancement/suanPan.sh DESTINATION . COMPONENT Support)
    set(DECOR ".")
    set(SUFFIX "")
    if (COMPILER_IDENTIFIER MATCHES "linux")
        set(DECOR ".so")
        set(SUFFIX "")
    elseif (COMPILER_IDENTIFIER MATCHES "mac")
        set(DECOR "")
        set(SUFFIX ".dylib")
        file(GLOB OPENBLAS_FILES Libs/${SP_EXTERNAL_LIB_PATH}/libopenblas*)
        install(PROGRAMS ${OPENBLAS_FILES} TYPE LIB COMPONENT Main)
    endif ()
    if (SP_BUILD_PARALLEL AND NOT SP_ENABLE_SYSLIB)
        file(GLOB TBB_FILES Libs/${SP_EXTERNAL_LIB_PATH}/libtbb*)
        install(PROGRAMS ${TBB_FILES} TYPE LIB COMPONENT Main)
    endif ()
    if (SP_ENABLE_MKL)
        if (SP_ENABLE_SHARED_MKL)
            file(GLOB MKL_FILES
                    ${MKL_ROOT}/lib/intel64/libmkl_core${DECOR}*${SUFFIX}
                    ${MKL_ROOT}/lib/intel64/libmkl_def${DECOR}*${SUFFIX}
                    ${MKL_ROOT}/lib/intel64/libmkl_avx*${DECOR}*${SUFFIX}
            )
            install(PROGRAMS ${MKL_FILES} TYPE LIB COMPONENT Main)
            if (SP_ENABLE_64BIT_INDEXING)
                file(GLOB MKL_FILES ${MKL_ROOT}/lib/intel64/libmkl_intel_ilp64${DECOR}*${SUFFIX})
            else ()
                file(GLOB MKL_FILES ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64${DECOR}*${SUFFIX})
            endif ()
            install(PROGRAMS ${MKL_FILES} TYPE LIB COMPONENT Main)
            if (SP_ENABLE_IOMP)
                file(GLOB MKL_FILES ${MKL_ROOT}/lib/intel64/libmkl_intel_thread${DECOR}*${SUFFIX})
            else ()
                file(GLOB MKL_FILES ${MKL_ROOT}/lib/intel64/libmkl_gnu_thread${DECOR}*${SUFFIX})
            endif ()
            install(PROGRAMS ${MKL_FILES} TYPE LIB COMPONENT Main)
        endif ()
        if (SP_ENABLE_IOMP)
            file(GLOB MKL_FILES ${IOMPPATH}/libiomp5${DECOR}*${SUFFIX})
            install(PROGRAMS ${MKL_FILES} TYPE LIB COMPONENT Main)
        endif ()
    endif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
        foreach (ONEAPI_PATH
                $ENV{ONEAPI_ROOT}/compiler/latest/lib #  since 2024
                $ENV{ONEAPI_ROOT}/compiler/latest/linux/compiler/lib/intel64_lin # prior to 2024
        )
            file(GLOB ONEAPI_FILES
                    ${ONEAPI_PATH}/libifcoremt${DECOR}*${SUFFIX}
                    ${ONEAPI_PATH}/libifport${DECOR}*${SUFFIX}
                    ${ONEAPI_PATH}/libimf${DECOR}*${SUFFIX}
                    ${ONEAPI_PATH}/libintlc${DECOR}*${SUFFIX}
                    ${ONEAPI_PATH}/libiomp5${DECOR}*${SUFFIX}
                    ${ONEAPI_PATH}/libsvml${DECOR}*${SUFFIX}
            )
            install(PROGRAMS ${ONEAPI_FILES} TYPE LIB COMPONENT Main)
        endforeach ()
    endif ()
    if (SP_ENABLE_MAGMA)
        file(GLOB MAGMA_FILES ${MAGMA_PATH}/lib/libmagma*.so)
        install(PROGRAMS ${MAGMA_FILES} TYPE LIB COMPONENT Main)
    endif ()
elseif (COMPILER_IDENTIFIER MATCHES "win")
    if (SP_BUILD_PARALLEL)
        file(GLOB TBB_FILES Libs/${SP_EXTERNAL_LIB_PATH}/libtbb*.dll)
        install(FILES ${TBB_FILES} DESTINATION bin COMPONENT Main)
    endif ()
    if (SP_ENABLE_MKL)
        file(GLOB MKL_FILES
                ${MKL_ROOT}/redist/intel64/mkl_rt*
                ${MKL_ROOT}/redist/intel64/mkl_core*
                ${MKL_ROOT}/redist/intel64/mkl_def*
                ${MKL_ROOT}/redist/intel64/mkl_avx*
                ${MKL_ROOT}/redist/intel64/mkl_intel_thread.*
                ${MKL_ROOT}/bin/mkl_rt*
                ${MKL_ROOT}/bin/mkl_core*
                ${MKL_ROOT}/bin/mkl_def*
                ${MKL_ROOT}/bin/mkl_avx*
                ${MKL_ROOT}/bin/mkl_intel_thread.*
        )
        install(FILES ${MKL_FILES} DESTINATION bin COMPONENT Main)
    endif ()
    file(GLOB DLL_FILES Libs/${SP_EXTERNAL_LIB_PATH}/lib*.dll)
    install(FILES ${DLL_FILES} DESTINATION bin COMPONENT Main)
elseif (COMPILER_IDENTIFIER MATCHES "vs")
    if (SP_BUILD_PARALLEL OR SP_ENABLE_MKL)
        file(GLOB TBB_FILES Libs/${SP_EXTERNAL_LIB_PATH}/tbb*.dll)
        install(FILES ${TBB_FILES} DESTINATION bin COMPONENT Main)
    endif ()
    if (SP_ENABLE_MKL)
        if (SP_ENABLE_SHARED_MKL)
            file(GLOB MKL_FILES
                    ${MKL_ROOT}/redist/intel64/mkl_core*
                    ${MKL_ROOT}/redist/intel64/mkl_def*
                    ${MKL_ROOT}/redist/intel64/mkl_avx*
                    ${MKL_ROOT}/redist/intel64/mkl_intel_thread.*
                    ${MKL_ROOT}/bin/mkl_core*
                    ${MKL_ROOT}/bin/mkl_def*
                    ${MKL_ROOT}/bin/mkl_avx*
                    ${MKL_ROOT}/bin/mkl_intel_thread.*
            )
            install(FILES ${MKL_FILES} DESTINATION bin COMPONENT Main)
        endif ()
        foreach (IOMP_DLL libifcoremd libiomp5md libmmd svml_dispmd)
            find_file(IOMP_${IOMP_DLL} ${IOMP_DLL}.dll PATHS ${MKL_ROOT}/../../ ${MKL_ROOT}/../../compiler/latest/bin/ REQUIRED NO_DEFAULT_PATH)
            install(FILES ${IOMP_${IOMP_DLL}} DESTINATION bin COMPONENT Main)
        endforeach ()
    else ()
        file(GLOB DLL_FILES Libs/win/lib*.dll)
        install(FILES ${DLL_FILES} DESTINATION bin COMPONENT Main)
    endif ()
    file(GLOB VC_FILES Libs/vc/*.dll)
    install(FILES ${VC_FILES} DESTINATION bin COMPONENT Main)
endif ()

message(STATUS "Flags and Dirs:")
message(STATUS "suanPan CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
if (CMAKE_Fortran_COMPILER)
    message(STATUS "suanPan Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
endif ()

message(STATUS "Link Dirs:")
get_property(SP_LINK_DIR DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY LINK_DIRECTORIES)
foreach (SDIR ${SP_LINK_DIR})
    message(STATUS "${SDIR}")
endforeach ()

message(STATUS "Include Dirs:")
get_property(SP_INCLUDE_DIR DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach (SDIR ${SP_INCLUDE_DIR})
    message(STATUS "${SDIR}")
endforeach ()

if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    file(READ "/etc/os-release" DISTRO_INFO)
    string(REGEX MATCH "fedora|ubuntu|debian" DIST ${DISTRO_INFO})

    if (DIST OR BUILD_PACKAGE)
        set(CPACK_PACKAGE_CONTACT "Theodore Chang")
        set(CPACK_PACKAGE_CHECKSUM "SHA256")
        set(CPACK_PACKAGE_ICON ${CMAKE_CURRENT_SOURCE_DIR}/Resource/suanPan-ua.svg)
        set(CPACK_PACKAGE_RELEASE 1)
        set(CPACK_PACKAGE_VENDOR "tlcfem")
        set(CPACK_PACKAGE_VERSION "4.0.0")
        set(CPACK_PACKAGE_DESCRIPTION "An Open Source, Parallel and Heterogeneous Finite Element Analysis Framework")
        set(CPACK_PACKAGE_HOMEPAGE_URL "https://github.com/TLCFEM/suanPan")
        set(CPACK_PACKAGE_NAME "${PROJECT_NAME}")
        set(CPACK_PACKAGE_FILE_NAME "${PROJECT_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_PACKAGE_RELEASE}.${CMAKE_SYSTEM_PROCESSOR}")
        set(CPACK_THREADS 0)
        set(CPACK_COMPONENTS_ALL Main)

        if ((DIST STREQUAL "ubuntu") OR (DIST STREQUAL "debian") OR (BUILD_PACKAGE MATCHES "DEB"))
            message(STATUS "Build DEB Package For Distribution: ${DIST}")
            set(CPACK_GENERATOR "DEB")
            set(CPACK_DEB_COMPONENT_INSTALL ON)
            set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
            set(CPACK_DEBIAN_COMPRESSION_TYPE "xz")
        elseif ((DIST STREQUAL "fedora") OR (BUILD_PACKAGE MATCHES "RPM"))
            message(STATUS "Build RPM Package For Distribution: ${DIST}")
            set(CPACK_GENERATOR "RPM")
            set(CPACK_RPM_COMPONENT_INSTALL ON)
            set(CPACK_RPM_PACKAGE_LICENSE "GPL-3.0")
        endif ()

        include(CPack)
    endif ()
endif ()

if (UNIX)
    get_filename_component(SP_BUILD_DIR "${CMAKE_BINARY_DIR}" NAME)
    if (SP_BUILD_DIR MATCHES "${CMAKE_BUILD_TYPE}")
        set(SP_SRC_JSON "${CMAKE_BINARY_DIR}/compile_commands.json")
        set(SP_DST_JSON "${CMAKE_BINARY_DIR}/../compile_commands.json")
        if (EXISTS "${SP_SRC_JSON}")
            execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink "${SP_SRC_JSON}" "${SP_DST_JSON}")
            message(STATUS "Created symlink: ${SP_DST_JSON} -> ${SP_SRC_JSON}")
        endif ()
    endif ()
endif ()