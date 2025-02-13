cmake_minimum_required(VERSION 3.14) 
# Set your compiler version
set(CMAKE_C_COMPILER gcc-13)
set(CMAKE_CXX_COMPILER g++-13)
project(GEM)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include(FetchContent)

# CMake options=================================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Fetch fmt ======================================================
FetchContent_Declare(fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG master
)
FetchContent_MakeAvailable(fmt)

# Fetch Armadillo=================================================================
FetchContent_Declare(
  armadillo
  GIT_REPOSITORY https://gitlab.com/conradsnicta/armadillo-code
  GIT_TAG 14.0.x
)

FetchContent_GetProperties(armadillo)
if(NOT armadillo_POPULATED)
  FetchContent_Populate(armadillo)
  add_library(armadillo INTERFACE)
  target_include_directories(armadillo INTERFACE ${armadillo_SOURCE_DIR}/include)
  # add_subdirectory(${armadillo_SOURCE_DIR} ${armadillo_BINARY_DIR})
endif()
# Fetch Eigen=====================================================
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
#Fetch Boost=================================================================
# FetchContent_Declare(
#  Boost
#  GIT_REPOSITORY https://github.com/boostorg/boost.git
#  GIT_TAG boost-1.85.0
#  GIT_SHALLOW TRUE
# )
# FetchContent_MakeAvailable(Boost)
# set(BOOST_INCLUDE_LIBRARIES thread filesystem system program_options asio algorithm iostreams math smart_ptr format)
# set(BOOST_ENABLE_CMAKE ON)

# #include(FetchContent)
# FetchContent_Declare(
#   Boost
#   GIT_REPOSITORY https://github.com/boostorg/boost.git
#   GIT_TAG boost-1.79.0
#   GIT_SHALLOW TRUE
# )
# FetchContent_MakeAvailable(Boost)

#find Boost=================================================================
set(BLA_VENDOR Intel10_64lp)
find_package(Boost 1.74.0 REQUIRED COMPONENTS program_options thread system filesystem)
include_directories(${Boost_INCLUDE_DIRS})
#find BLAS/LAPACK=================================================================
# find_package(BLAS REQUIRED)
# find_package(LAPACK REQUIRED)
#Set MKL=================================================================
set(MKLROOT "/opt/intel/compilers_and_libraries_2020.2.254/linux/mkl")
message("MKLROOT: ${MKLROOT}")
message(STATUS "Setting BLAS vendor to Intel10_64lp.")

# Specify include directories for MKL
include_directories(${MKLROOT}/include)
message(STATUS "Setting include directory for MKL to ${MKLROOT}/include.")

# Specify library directories for MKL
link_directories(${MKLROOT}/lib/intel64)
message(STATUS "Setting library directory for MKL to ${MKLROOT}/lib/intel64.")

# Define EIGEN_USE_MKL_ALL for Eigen library
add_definitions(-DEIGEN_USE_MKL_ALL)
# Define ARMA to use MKL
add_definitions(-DARMA_USE_LAPACK -DARMA_USE_BLAS -DARMA_DONT_USE_WRAPPER)

message(STATUS "Defining EIGEN_USE_MKL_ALL for Eigen  and Armadillo libraries.")

# Set MKL libraries
set(MKL_LIBRARIES 
    mkl_intel_lp64       
    mkl_sequential    
    mkl_core            
    pthread
    m
    dl
)


#SuiteSparse=================================================================
#If SuiteSparse is installed in default directory uncomment below line and comment all the lines related to the SuiteSparse
#find_package(SuiteSparse REQUIRED)
# Define SuiteSparse installation directory
# Include the SuiteSparse headers compiled with suitesparse_7.6.0_no_omp replace your path
include_directories("/home/ssalehinasab/local/suitesparse_7.6.0_no_omp/include/suitesparse/")

# Link the SuiteSparse libraries
link_directories("/home/ssalehinasab/local/suitesparse_7.6.0_no_omp/lib64")

# Add the SuiteSparse libraries to the target
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

#=================================================================
file(GLOB SRC "src/*.cpp")
file(GLOB CC_SOURCES "thirdparty/plink-2.0/*.cc")
file(GLOB C_SOURCES 
    "thirdparty/zstd-1.5.5/lib/common/*.c"
    "thirdparty/zstd-1.5.5/lib/compress/*.c"
    "thirdparty/zstd-1.5.5/lib/decompress/*.c"
    "thirdparty/libdeflate-1.18/lib/*.c"
    "thirdparty/libdeflate-1.18/lib/x86/*.c"
    "thirdparty/libdeflate-1.18/lib/arm/*.c"
)
# Concatenate the two lists of source files
set(SOURCES ${SRC} ${CC_SOURCES} ${C_SOURCES})
#=================================================================
# Add the executable with the concatenated source files
add_executable(GEM ${SOURCES} )
target_link_libraries(GEM
  PRIVATE
    fmt::fmt
    # BLAS::BLAS 
    # LAPACK::LAPACK
    eigen
    armadillo 
    Boost::thread
    Boost::system
    Boost::filesystem
    Boost::program_options
    ${SUITESPARSE_LIBRARIES}
    ${MKL_LIBRARIES}
)

# Set compiler optimization flags 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -funroll-loops -ffast-math  -fno-math-errno -DNDEBUG") #-Wall -DMKL_ILP64 -march=native -m64 -fopenmp -mfma
