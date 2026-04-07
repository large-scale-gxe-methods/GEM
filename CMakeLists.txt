cmake_minimum_required(VERSION 3.14)
# Set the compiler
set(CMAKE_C_COMPILER gcc-13)
set(CMAKE_CXX_COMPILER g++-13)
project(GEM LANGUAGES C CXX)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
# =========================================================
# Options
# =========================================================
option(BUILD_STATIC "Build static executable" OFF)
# option(USE_MKL "Use Intel MKL instead of BLAS/LAPACK" OFF)
set(BLAS_BACKEND "MKL" CACHE STRING "BLAS backend: AUTO, MKL, OPENBLAS")
set_property(CACHE BLAS_BACKEND PROPERTY STRINGS AUTO MKL OPENBLAS)
# =========================================================
# C++ standard
# =========================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# =========================================================
# Static vs Dynamic
# =========================================================
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g")
if(BUILD_STATIC)
    message(STATUS "Building STATIC executable")
    set(BUILD_SHARED_LIBS OFF)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -no-pie -static")
else()
    message(STATUS "Building DYNAMIC executable")
endif()
# =========================================================
# fmt
# =========================================================
include(FetchContent)

FetchContent_Declare(fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG 11.0.2
)
FetchContent_MakeAvailable(fmt)
# =========================================================
# spdlog
# =========================================================
set(SPDLOG_BUILD_SHARED OFF CACHE BOOL "" FORCE)
FetchContent_Declare(spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG v1.15.0
)
FetchContent_MakeAvailable(spdlog)
# =========================================================
# Fetch Eigen
# =========================================================
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen
    GIT_TAG master
)
FetchContent_GetProperties(eigen)
if(NOT eigen_POPULATED)
    FetchContent_Populate(eigen)
    add_library(eigen INTERFACE)
    target_include_directories(eigen INTERFACE ${eigen_SOURCE_DIR})
endif()
message(STATUS "Eigen source dir = ${eigen_SOURCE_DIR}")
# =========================================================
# Armadillo
# =========================================================
FetchContent_Declare(
  armadillo
  GIT_REPOSITORY https://gitlab.com/conradsnicta/armadillo-code
  GIT_TAG 14.0.x
)
FetchContent_GetProperties(armadillo)
if(NOT armadillo_POPULATED)
  FetchContent_Populate(armadillo)
  include_directories(${armadillo_SOURCE_DIR}/include)
  # add_subdirectory(${armadillo_SOURCE_DIR} ${armadillo_BINARY_DIR})
endif()
# =========================================================
# Boost (required)
# =========================================================
if(BUILD_STATIC)
    set(Boost_USE_STATIC_LIBS ON)
else()
    set(Boost_USE_STATIC_LIBS OFF)
endif()

find_package(Boost REQUIRED COMPONENTS
    program_options
    thread
    system
    filesystem
)
# =========================================================
# BLAS Backend Selection
# =========================================================
if(BLAS_BACKEND STREQUAL "MKL")

    message(STATUS "Configuring MKL backend")

    if(DEFINED ENV{MKLROOT})
        set(MKLROOT $ENV{MKLROOT})
        message(STATUS "Using exported MKL path: ${MKLROOT}")
    elseif()
        set(MKLROOT "/opt/intel/compilers_and_libraries_2020.2.254/linux/mkl")
        message(STATUS "Using default MKL path: ${MKLROOT}")
    else()
        message(FATAL_ERROR
            "MKL requested but MKLROOT not found.\n"
            "Please set MKLROOT, e.g.:\n"
            "export MKLROOT=/opt/intel/oneapi/mkl/latest"
        )
    endif()

    if(BUILD_STATIC)
        message(STATUS "Using STATIC MKL")
        set(BLAS_LIBS
            -Wl,--start-group
            ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a
            ${MKLROOT}/lib/intel64/libmkl_sequential.a
            ${MKLROOT}/lib/intel64/libmkl_core.a
            -Wl,--end-group
            -Wl,-Bstatic -ldl -lpthread -lm
        )
    else()
        set(BLAS_LIBS
            mkl_intel_lp64
            mkl_sequential
            mkl_core
            pthread m dl
        )
    endif()

elseif(BLAS_BACKEND STREQUAL "OPENBLAS")

    message(STATUS "Configuring OpenBLAS backend")

    if(BUILD_STATIC)
        message(STATUS "BUILD_STATIC is ON → trying static OpenBLAS")

        find_library(OPENBLAS_STATIC
            NAMES openblas
            PATH_SUFFIXES lib lib64
        )

        if(OPENBLAS_STATIC AND OPENBLAS_STATIC MATCHES "\\.a$")
            message(STATUS "Found static OpenBLAS: ${OPENBLAS_STATIC}")

           set(BLAS_LIBS
            # -Wl,-Bstatic
            ${OPENBLAS_STATIC}
            # -Wl,-Bdynamic #dynamic libs
            # gfortran
            # dl
            # pthread
            # m
            -Wl,-Bstatic -lgfortran -ldl -lpthread -lm
        )

        else()
            message(STATUS "System OpenBLAS not found. Using local OpenBLAS.")

            set(OPENBLAS_LOCAL "$ENV{HOME}/local/openblas")

            set(BLAS_LIBS ${OPENBLAS_LOCAL}/lib/libopenblas.a)
            include_directories(${OPENBLAS_LOCAL}/include)

            if (NOT EXISTS ${BLAS_LIBS})
                message(FATAL_ERROR "Local OpenBLAS not found in ${OPENBLAS_LOCAL}")
            endif()

        endif()

    else()

        message(STATUS "BUILD_STATIC is OFF → using shared OpenBLAS")

        find_library(OPENBLAS_LIB
            NAMES openblas
            PATHS /usr/lib64
        )

        if(NOT OPENBLAS_LIB)
            message(FATAL_ERROR "System OpenBLAS not found")
        endif()

        message(STATUS "Using system OpenBLAS: ${OPENBLAS_LIB}")
        set(BLAS_LIBS ${OPENBLAS_LIB})

    endif()

elseif(BLAS_BACKEND STREQUAL "AUTO")

    message(STATUS "Configuring AUTO BLAS detection")

    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    set(BLAS_LIBS ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

else()
    message(FATAL_ERROR "Invalid BLAS_BACKEND: ${BLAS_BACKEND}")
endif()

# =========================================================
# SuiteSparse (required)
# =========================================================
# Uncomment the line below and comment out the following lines if SuiteSparse is installed in the system root directory
# find_package(SuiteSparse REQUIRED)
# If SuiteSparse is not installed in the system root directory, set the paths manually
include_directories("$ENV{HOME}/local/suitesparse_7.6.0_no_omp/include/suitesparse")
link_directories("$ENV{HOME}/local/suitesparse_7.6.0_no_omp/lib64")
set(SUITESPARSE_LIBRARIES
    libcholmod.a 
    libcxsparse.a
    libspqr.a 
    libumfpack.a 
    libcamd.a 
    libccolamd.a 
    libcolamd.a 
    libamd.a 
    libsuitesparseconfig.a 
)

# =========================================================
# Sources
# =========================================================
file(GLOB CXX_SOURCES "src/*.cpp")

file(GLOB C_SOURCES 
    "thirdparty/zstd-1.5.5/lib/common/*.c"
    "thirdparty/zstd-1.5.5/lib/compress/*.c"
    "thirdparty/zstd-1.5.5/lib/decompress/*.c"
    "thirdparty/libdeflate-1.18/lib/*.c"
    "thirdparty/libdeflate-1.18/lib/x86/*.c"
    "thirdparty/libdeflate-1.18/lib/arm/*.c"
)

file(GLOB CC_SOURCES "thirdparty/plink-2.0/*.cc")

# =========================================================
# Executable
# =========================================================
add_executable(GEM2
    ${CXX_SOURCES}
    ${C_SOURCES}
    ${CC_SOURCES}
)
# Armadillo always uses BLAS/LAPACK
target_compile_definitions(GEM2 PRIVATE
    ARMA_USE_BLAS
    ARMA_USE_LAPACK
    ARMA_DONT_USE_WRAPPER
)
if(BLAS_BACKEND STREQUAL "MKL")
    target_compile_definitions(GEM2 PRIVATE EIGEN_USE_MKL_ALL)

elseif(BLAS_BACKEND STREQUAL "OPENBLAS")
    target_compile_definitions(GEM2 PRIVATE
        EIGEN_USE_BLAS
        EIGEN_USE_LAPACKE
    )
endif()

# =========================================================
# Link libraries
# =========================================================
target_link_libraries(GEM2 PRIVATE
    eigen
    fmt::fmt
    spdlog::spdlog
    Boost::program_options
    Boost::thread
    Boost::system
    Boost::filesystem
    ${SUITESPARSE_LIBRARIES}
    ${BLAS_LIBS}
)



