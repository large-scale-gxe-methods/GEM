CC=gcc  # Do not change this for libdeflate package.
CXX=g++
CXXFLAGS = -std=c++11 -O2 -g

BOOST_INC=/usr/local/include 
Eigen_INC=/home/cpan/eigen/eigen-3.4.0
BOOST_LIB_TD=/usr/local/lib/libboost_thread.a
BOOST_LIB_PO=/usr/local/lib/libboost_program_options.a
BOOST_LIB_SY=/usr/local/lib/libboost_system.a
BOOST_LIB_FS=/usr/local/lib/libboost_filesystem.a

DYN_LIBS=-lmkl_gf_lp64 -lmkl_sequential -lmkl_core -ldl -lpthread -lm -larmadillo


# Uncomment line below for LAPACK/BLAS only and comment other lines beginning with LIBS.
#LIBS= -llapack -lblas -ldl -pthread -lm -lboost_program_options -lboost_thread-mt -lboost_system -lboost_filesystem

#LIBS= -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -ldl -lpthread -lm -lboost_program_options -lboost_thread -lboost_system -lboost_filesystem

C_SOURCES= $(wildcard ../thirdparty/zstd-1.4.5/lib/common/*.c) $(wildcard ../thirdparty/zstd-1.4.5/lib/compress/*.c) $(wildcard ../thirdparty/zstd-1.4.5/lib/decompress/*.c) $(wildcard ../thirdparty/zstd-1.4.5/lib/dictBuilder/*.c) $(wildcard ../thirdparty/libdeflate-1.6/lib/*.c) $(wildcard ../thirdparty/libdeflate-1.6/lib/x86/*.c) $(wildcard ../thirdparty/libdeflate-1.6/lib/arm/*.c)
CXX_SOURCES= MatrixUtils.cpp ReadParameters.cpp TimeUtils.cpp BinaryEUtils.cpp ReadBGEN.cpp ReadPGEN.cpp ReadBed.cpp GEM.cpp $(wildcard ../thirdparty/plink-2.0/*.cpp)
CXX_OBJECTS= $(CXX_SOURCES:.cpp=.o)
C_OBJECTS= $(C_SOURCES:.c=.o)
OBJS= $(CXX_OBJECTS) $(C_OBJECTS)

#targets
all: GEM

GEM: $(CXX_OBJECTS) $(C_OBJECTS)
	$(CXX) $^ $(CXXFLAGS) $(BOOST_LIB_TD) $(BOOST_LIB_PO) $(BOOST_LIB_SY) $(BOOST_LIB_FS) -o $@ $(DYN_LIBS)

%.o: %.c
	$(CC) -O2 -g -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -I$(BOOST_INC) -I$(Eigen_INC)

.PHONY: clean
clean:
	rm -f $(OBJS)
