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

include_directories(${ROOT})
include_directories(Include)
include_directories(Include/metis)
include_directories(Include/fmt/include)

message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

set(BUILD_PACKAGE "" CACHE STRING "DEB OR RPM")

option(BUILD_SHARED_LIBS "Build shared libraries." OFF)
option(BUILD_DLL_EXAMPLE "Build dynamic linked library examples." OFF)
option(BUILD_MULTITHREAD "Build with multi-threaded support via TBB." OFF)
option(USE_SUPERLUMT "Use multi-threaded SuperLU. Note the performance may not be better than the sequential version." OFF)
option(USE_VTK "Enable visualisation via VTK. Note external VTK libraries need to be compiled in advance." OFF)
option(USE_HDF5 "Enable recording results in HDF5 format." ON)
option(USE_AVX "Enable AVX support." OFF)
option(USE_AVX2 "Enable AVX2 support." ON)
option(USE_AVX512 "Enable AVX512 support." OFF)
option(USE_MKL "Use Intel MKL instead of OpenBLAS." OFF)
if (LINUX AND CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    option(USE_AOCL "Use AMD Optimizing CPU Libraries instead of OpenBLAS." OFF)
endif ()
option(USE_MIMALLOC "Use mimalloc instead of default memory allocator." OFF)
option(USE_SYS_LIB "Use libraries installed on the system instead of the bundled ones." OFF)
if (USE_MKL)
    option(USE_INTEL_OPENMP "Use Intel OpenMP implementation on Linux and macOS" ON)
    option(LINK_DYNAMIC_MKL "Link dynamic Intel MKL libraries." ON)
    option(USE_MPI "Enable MPI based global solvers." OFF)
else ()
    set(USE_MPI OFF CACHE BOOL "" FORCE)
    set(CUSTOM_OPENBLAS "" CACHE PATH "The path that contains the custom OpenBLAS library. If not set, the bundled OpenBLAS will be used.")
endif ()
if (USE_AOCL)
    set(AOCL_BLIS_PATH "${CMAKE_SOURCE_DIR}/Libs/linux/libblis-mt.a" CACHE FILEPATH "AOCL libblis path.")
    set(AOCL_FLAME_PATH "${CMAKE_SOURCE_DIR}/Libs/linux/libflame.a" CACHE FILEPATH "AOCL libflame path.")
    set(AOCL_UTILS_PATH "${CMAKE_SOURCE_DIR}/Libs/linux/libaoclutils.a" CACHE FILEPATH "AOCL libaoclutils path.")
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
        option(USE_CUDA "Enable GPU based global solvers via CUDA." OFF)
    endif ()
elseif (CMAKE_SYSTEM_NAME MATCHES "Linux") # LINUX PLATFORM
    set(SP_EXTERNAL_LIB_PATH "linux")
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU") # GNU GCC COMPILER
        set(COMPILER_IDENTIFIER "gcc-linux")
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM") # Intel COMPILER icpx
        set(COMPILER_IDENTIFIER "clang-linux")
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel") # Intel COMPILER Classic icc
        set(COMPILER_IDENTIFIER "gcc-linux")
        message(STATUS "Classic Intel compiler icc has incomplete CPP20 support, if it fails to compile please use another compiler.")
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang") # Clang COMPILER
        set(COMPILER_IDENTIFIER "clang-linux")
    endif ()
    option(USE_CUDA "Enable GPU based global solvers via CUDA." OFF)
elseif (CMAKE_SYSTEM_NAME MATCHES "Darwin") # MAC PLATFORM
    set(SP_EXTERNAL_LIB_PATH "mac")
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU") # GNU GCC COMPILER
        set(COMPILER_IDENTIFIER "gcc-mac")
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(COMPILER_IDENTIFIER "clang-mac")
        link_directories(/usr/local/opt/libomp/lib/)
        message(STATUS "On macOS, make sure llvm and libomp are installed.")
        message(STATUS "brew install llvm libomp")
    endif ()
endif ()

if (COMPILER_IDENTIFIER MATCHES "unknown")
    message(FATAL_ERROR "Cannot identify the compiler available, please use GCC or MSVC or Intel.")
endif ()

if (NOT CUSTOM_OPENBLAS STREQUAL "")
    link_directories(${CUSTOM_OPENBLAS})
elseif (USE_SYS_LIB)
    find_package(OpenBLAS REQUIRED)
endif ()

if (NOT USE_SYS_LIB)
    link_directories(Libs/${SP_EXTERNAL_LIB_PATH})
endif ()

if (USE_SUPERLUMT)
    message(WARNING "Current SuperLU MT library may contain bugs. Disabling it.")
    set(USE_SUPERLUMT OFF CACHE BOOL "" FORCE)
    # add_compile_definitions(SUANPAN_SUPERLUMT)
endif ()

if (USE_MPI)
    set(MPI_EXECUTABLE_SUFFIX "")
    find_package(MPI REQUIRED)
    add_compile_definitions(SUANPAN_MPI MPICH_SKIP_MPICXX)
    include_directories(${MPI_C_INCLUDE_DIRS})
    string(TOLOWER "${MPI_C_INCLUDE_DIRS}" MPI_C_INCLUDE_DIRS_LOWER)
    if (MPI_C_INCLUDE_DIRS_LOWER MATCHES "openmpi")
        set(MKL_MPI openmpi)
    elseif (MPI_C_INCLUDE_DIRS_LOWER MATCHES "mpich")
        set(MKL_MPI mpich)
    endif ()
endif ()

if (USE_MKL)
    if (USE_CUDA)
        option(USE_MAGMA "Enable GPU based global solvers via MAGMA." OFF)
    endif ()
    set(MKLROOT "" CACHE PATH "MKL library path which contains /include and /lib folders.")
    if (COMPILER_IDENTIFIER MATCHES "IntelLLVM")
        set(USE_INTEL_OPENMP ON CACHE BOOL "" FORCE)
    endif ()
    if (NOT LINK_DYNAMIC_MKL)
        set(MKL_LINK static)
    endif ()
    if (NOT USE_INTEL_OPENMP)
        set(MKL_THREADING gnu_thread)
    endif ()
    set(MKL_INTERFACE lp64)
    find_package(MKL REQUIRED HINTS ${MKLROOT} "/opt/intel/oneapi/mkl/latest" "C:/Program Files (x86)/Intel/oneAPI/mkl/latest")
    add_compile_definitions(SUANPAN_MKL)
    # add_compile_definitions(ARMA_USE_MKL_ALLOC)
    message(STATUS "MKL_ROOT: ${MKL_ROOT}")
elseif (FORTRAN_STATUS AND CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    message(WARNING "Since Intel compilers are used, why not enabling MKL?")
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
    add_compile_definitions(SUANPAN_AOCL)
endif ()

if (USE_CUDA)
    if (POLICY CMP0146)
        cmake_policy(SET CMP0146 OLD)
    endif ()
    find_package(CUDA PATHS ${CUDA_PATH})
    if (NOT CUDA_FOUND)
        set(CUDA_PATH "" CACHE PATH "CUDA library path which contains /include folder.")
        message(FATAL_ERROR "CUDA library is not found, please indicate its path.")
    endif ()
    add_compile_definitions(SUANPAN_CUDA)
    include_directories(${CUDA_INCLUDE_DIRS})
    link_libraries(${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_cusolver_LIBRARY} ${CUDA_cusparse_LIBRARY})
    if (CUDA_VERSION VERSION_GREATER_EQUAL "12.3" AND COMPILER_IDENTIFIER MATCHES "vs")
        add_compile_options(/wd4996)
    endif ()
endif ()

if (USE_MAGMA)
    set(MAGMAROOT "" CACHE PATH "Magma library path which contains /include and /lib folders.")
    find_file(MAGMA_HEADER NAMES magma.h PATHS ${MAGMAROOT}/include)
    if (MAGMA_HEADER MATCHES "MAGMA_HEADER-NOTFOUND")
        message(FATAL_ERROR "The <magma.h> is not found under the path: ${MAGMAROOT}/include.")
    endif ()
    include_directories(${MAGMAROOT}/include)
    link_directories(${MAGMAROOT}/lib)
    link_libraries(magma magma_sparse)
    add_compile_definitions(SUANPAN_MAGMA)
endif ()

if (USE_VTK)
    find_package(VTK PATHS ${VTK_PATH})
    if (NOT VTK_FOUND)
        set(VTK_PATH "" CACHE PATH "VTK library path which contains /include folder.")
        message(FATAL_ERROR "VTK library is not found, please indicate its path.")
    endif ()
    add_compile_definitions(SUANPAN_VTK)
endif ()

if (USE_HDF5)
    add_compile_definitions(SUANPAN_HDF5)
    if (USE_VTK)
        string(REGEX REPLACE "/lib6?4?/cmake/vtk" "/include/vtk" VTK_INCLUDE ${VTK_DIR}) # on linux
        string(REGEX REPLACE "\\\\lib6?4?\\\\cmake\\\\vtk" "\\\\include\\\\vtk" VTK_INCLUDE ${VTK_INCLUDE}) # on windows
        include_directories(${VTK_INCLUDE}/vtkhdf5 ${VTK_INCLUDE}/vtkhdf5/src ${VTK_INCLUDE}/vtkhdf5/hl/src)
        find_file(HDF5_HEADER NAMES hdf5.h PATHS ${VTK_INCLUDE}/vtkhdf5/src)
        if (HDF5_HEADER MATCHES "HDF5_HEADER-NOTFOUND")
            message(FATAL_ERROR "The <hdf5.h> is not found in the include directories.")
        else ()
            message(STATUS "Found HDF5 header: ${HDF5_HEADER}")
        endif ()
    elseif (USE_SYS_LIB)
        find_package(HDF5 REQUIRED COMPONENTS C HL)
        include_directories(${HDF5_INCLUDE_DIRS})
        link_libraries(${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
    else ()
        include_directories(Include/hdf5 Include/hdf5-${SP_EXTERNAL_LIB_PATH})
        if (COMPILER_IDENTIFIER MATCHES "vs")
            link_libraries(libhdf5_hl libhdf5 shlwapi)
        else ()
            link_libraries(hdf5_hl hdf5)
        endif ()
    endif ()
else ()
    add_compile_definitions(ARMA_DONT_USE_HDF5)
endif ()

if (BUILD_MULTITHREAD)
    message(STATUS "USING TBB LIBRARY")
    add_compile_definitions(SUANPAN_MT)
    if (USE_SYS_LIB)
        find_package(TBB REQUIRED)
        include_directories(${TBB_INCLUDE_DIRS})
        link_libraries(TBB::tbb TBB::tbbmalloc TBB::tbbmalloc_proxy)
    else ()
        if (COMPILER_IDENTIFIER MATCHES "gcc-win")
            link_libraries(tbb12)
        else ()
            link_libraries(tbb)
        endif ()
        option(USE_TBB_ALLOC "Use tbb memory allocator. Enable if no other allocators will be used." OFF)
        if (USE_TBB_ALLOC)
            # for armadillo to use tbb allocator
            add_compile_definitions(ARMA_USE_TBB_ALLOC)
            include_directories(Include/oneapi) # because armadillo assumes oneapi be in the include path
            link_libraries(tbbmalloc tbbmalloc_proxy)
        endif ()
    endif ()
endif ()

if (BUILD_SHARED_LIBS)
    message(STATUS "BUILD SHARED LIBRARY")
else ()
    message(STATUS "BUILD STATIC LIBRARY")
endif ()

if (USE_AVX512)
    add_compile_definitions(SUANPAN_AVX512)
elseif (USE_AVX2)
    add_compile_definitions(SUANPAN_AVX2)
elseif (USE_AVX)
    add_compile_definitions(SUANPAN_AVX)
endif ()

if (COMPILER_IDENTIFIER MATCHES "vs")
    unset(TEST_COVERAGE CACHE)

    link_directories(Libs/win)

    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MP /openmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /openmp /EHsc")

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /nowarn /MP /Qopenmp /Qparallel /fpp /names:lowercase /assume:underscore /libs:dll /threads")

    if (USE_AVX512)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /arch:AVX512")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX512")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /arch:AVX2")
    elseif (USE_AVX2)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /arch:AVX2")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /arch:AVX2")
    elseif (USE_AVX)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /arch:AVX")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /arch:AVX")
    endif ()
else ()
    if (BUILD_SHARED_LIBS)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")
    endif ()

    link_libraries(dl pthread gfortran)

    find_library(HAS_QUADMATH quadmath)
    if (HAS_QUADMATH)
        link_libraries(quadmath)
    endif ()

    if (CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
        link_libraries(stdc++)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffp-model=precise -fexceptions -fiopenmp")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffp-model=precise -fexceptions -fiopenmp")
    elseif (COMPILER_IDENTIFIER MATCHES "clang-mac")
        include_directories("/usr/local/include" "/usr/local/opt/llvm/include")
        link_directories("/usr/local/lib" "/usr/local/opt/llvm/lib")
        link_libraries(omp)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fexceptions -Xpreprocessor -fopenmp")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fexceptions -Xpreprocessor -fopenmp")
    else ()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fexceptions -fopenmp")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fexceptions -fopenmp")
    endif ()

    if (USE_AVX512)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512f")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx512f")
    elseif (USE_AVX2)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx2")
    elseif (USE_AVX)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx")
    endif ()

    option(TEST_COVERAGE "TEST CODE COVERAGE USING GCOV" OFF)

    if (TEST_COVERAGE AND COMPILER_IDENTIFIER MATCHES "gcc") # only report coverage with gcc
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
        link_libraries(gcov)
    endif ()

    if (CMAKE_BUILD_TYPE MATCHES "Debug")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
        option(USE_ASAN "USE ADDRESS SANITIZER" OFF)
        if (USE_ASAN)
            message(STATUS "Using the address sanitizer with flags: -fsanitize=address,leak,undefined")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address,leak,undefined")
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address,leak,undefined")
        endif ()
    endif ()

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w -fallow-argument-mismatch")
    if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp -qopenmp")

        if (USE_AVX512)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -arch AVX2")
        elseif (USE_AVX2)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -arch AVX2")
        elseif (USE_AVX)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -arch AVX")
        endif ()
    else ()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fopenmp")

        if (USE_AVX512)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mavx512f")
        elseif (USE_AVX2)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mavx2")
        elseif (USE_AVX)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mavx")
        endif ()
    endif ()
endif ()

if (USE_MKL)
    if (COMPILER_IDENTIFIER MATCHES "win")
        link_libraries(mkl_rt)
    else ()
        link_libraries(MKL::MKL)
    endif ()
elseif (COMPILER_IDENTIFIER MATCHES "vs")
    link_libraries(libopenblas)
elseif (USE_AOCL)
    link_libraries(${AOCL_FLAME_PATH} ${AOCL_BLIS_PATH} ${AOCL_UTILS_PATH})
else ()
    link_libraries(openblas)
    if (COMPILER_IDENTIFIER MATCHES "linux")
        set(CMAKE_EXE_LINKER_FLAGS "-Wl,-init,gotoblas_init")
    endif ()
endif ()

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

set_property(TARGET ${PROJECT_NAME} PROPERTY ENABLE_EXPORTS 1)

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

add_subdirectory(Toolbox/lis-src)
add_subdirectory(Toolbox/metis-src)
target_link_libraries(${PROJECT_NAME} lis metis)

if (USE_SUPERLUMT)
    add_subdirectory(Toolbox/superlumt-src)
    target_link_libraries(${PROJECT_NAME} superlumt)
else ()
    add_subdirectory(Toolbox/superlu-src)
    target_link_libraries(${PROJECT_NAME} superlu)
endif ()

if (USE_MPI)
    add_subdirectory(MPI)
    target_link_libraries(${PROJECT_NAME} MPI::MPI_C MPI::MPI_CXX MPI::MPI_Fortran)
endif ()

if (USE_MIMALLOC)
    message(STATUS "USING MIMALLOC LIBRARY")
    include(FetchContent)
    FetchContent_Declare(mimalloc
            GIT_REPOSITORY https://github.com/microsoft/mimalloc
            GIT_TAG v2.1.9)
    FetchContent_MakeAvailable(mimalloc)
    add_dependencies(${PROJECT_NAME} mimalloc)
endif ()

if (USE_VTK)
    target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})
    if (VTK_VERSION VERSION_LESS "8.90.0")
        include(${VTK_USE_FILE})
    else ()
        vtk_module_autoinit(TARGETS ${PROJECT_NAME} MODULES ${VTK_LIBRARIES})
    endif ()
endif ()

if (BUILD_DLL_EXAMPLE)
    add_subdirectory(Developer/Element)
    add_subdirectory(Developer/Material)
    add_subdirectory(Developer/Modifier)
    add_subdirectory(Developer/ModuleBundle)
    add_subdirectory(Developer/Section)
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
    if (USE_MKL AND MKL_ROOT MATCHES "(oneapi|oneAPI)")
        if (LINK_DYNAMIC_MKL)
            file(GLOB MKL_FILES ${MKL_ROOT}/lib/intel64/libmkl_core${DECOR}*${SUFFIX})
            install(PROGRAMS ${MKL_FILES} DESTINATION lib)
            file(GLOB MKL_FILES ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64${DECOR}*${SUFFIX})
            install(PROGRAMS ${MKL_FILES} DESTINATION lib)
            file(GLOB MKL_FILES ${MKL_ROOT}/lib/intel64/libmkl_def${DECOR}*${SUFFIX})
            install(PROGRAMS ${MKL_FILES} DESTINATION lib)
            file(GLOB MKL_FILES ${MKL_ROOT}/lib/intel64/libmkl_avx*${DECOR}*${SUFFIX})
            install(PROGRAMS ${MKL_FILES} DESTINATION lib)
            if (USE_INTEL_OPENMP)
                file(GLOB MKL_FILES ${MKL_ROOT}/lib/intel64/libmkl_intel_thread${DECOR}*${SUFFIX})
                install(PROGRAMS ${MKL_FILES} DESTINATION lib)
            else ()
                file(GLOB MKL_FILES ${MKL_ROOT}/lib/intel64/libmkl_gnu_thread${DECOR}*${SUFFIX})
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
    if (USE_MKL AND MKL_ROOT MATCHES "(oneapi|oneAPI)")
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
    if (USE_MKL AND MKL_ROOT MATCHES "(oneapi|oneAPI)")
        if (LINK_DYNAMIC_MKL)
            file(GLOB MKL_FILES
                    ${MKL_ROOT}/redist/intel64/mkl_core*
                    ${MKL_ROOT}/redist/intel64/mkl_def*
                    ${MKL_ROOT}/redist/intel64/mkl_avx*
                    ${MKL_ROOT}/bin/mkl_core*
                    ${MKL_ROOT}/bin/mkl_def*
                    ${MKL_ROOT}/bin/mkl_avx*
            )
            install(FILES ${MKL_FILES} DESTINATION bin)
            file(GLOB MKL_FILES
                    ${MKL_ROOT}/redist/intel64/mkl_intel_thread.*
                    ${MKL_ROOT}/bin/mkl_intel_thread.*
            )
            install(FILES ${MKL_FILES} DESTINATION bin)
        endif ()
        if (USE_INTEL_OPENMP)
            find_file(libifcoremd libifcoremd.dll PATHS ${MKL_ROOT}/../../ ${MKL_ROOT}/../../compiler/latest/bin/ REQUIRED)
            install(FILES ${libifcoremd} DESTINATION bin)
            find_file(libiomp5md libiomp5md.dll PATHS ${MKL_ROOT}/../../ ${MKL_ROOT}/../../compiler/latest/bin/ REQUIRED)
            install(FILES ${libiomp5md} DESTINATION bin)
            find_file(libmmd libmmd.dll PATHS ${MKL_ROOT}/../../ ${MKL_ROOT}/../../compiler/latest/bin/ REQUIRED)
            install(FILES ${libmmd} DESTINATION bin)
            find_file(svml_dispmd svml_dispmd.dll PATHS ${MKL_ROOT}/../../ ${MKL_ROOT}/../../compiler/latest/bin/ REQUIRED)
            install(FILES ${svml_dispmd} DESTINATION bin)
        endif ()
    else ()
        file(GLOB DLL_FILES Libs/win/lib*.dll)
        install(FILES ${DLL_FILES} DESTINATION bin)
    endif ()
    file(GLOB VC_FILES Libs/vc/*.dll)
    install(FILES ${VC_FILES} DESTINATION bin)
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
