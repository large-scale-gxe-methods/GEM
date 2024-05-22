CC=gcc  # Do not change this for libdeflate package.
CXX=g++
CXXFLAGS = -std=c++11 -O2 -g

BOOST_INC=../boost 
Eigen_INC=../eigen
BOOST_LIB=-lboost_thread -lboost_program_options -lboost_system -lboost_filesystem

DYN_LIBS=-lmkl_gf_lp64 -lmkl_sequential -lmkl_core -ldl -lpthread -lm


# Uncomment line below for LAPACK/BLAS only and comment other lines beginning with LIBS.
#LIBS= -llapack -lblas -ldl -pthread -lm -lboost_program_options -lboost_thread-mt -lboost_system -lboost_filesystem

#LIBS= -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -ldl -lpthread -lm -lboost_program_options -lboost_thread -lboost_system -lboost_filesystem

C_SOURCES= $(wildcard ../thirdparty/zstd-1.5.5/lib/common/*.c) $(wildcard ../thirdparty/zstd-1.5.5/lib/compress/*.c) $(wildcard ../thirdparty/zstd-1.5.5/lib/decompress/*.c) $(wildcard ../thirdparty/libdeflate-1.18/lib/*.c) $(wildcard ../thirdparty/libdeflate-1.18/lib/x86/*.c) $(wildcard ../thirdparty/libdeflate-1.18/lib/arm/*.c)
CXX_SOURCES= MatrixUtils.cpp ReadParameters.cpp TimeUtils.cpp BinaryEUtils.cpp ReadBGEN.cpp ReadPGEN.cpp ReadBed.cpp GEM.cpp 
CC_SOURCES=$(wildcard ../thirdparty/plink-2.0/*.cc)
CXX_OBJECTS= $(CXX_SOURCES:.cpp=.o)
C_OBJECTS= $(C_SOURCES:.c=.o)
CC_OBJ= $(CC_SOURCES:.cc=.o)
OBJS= $(CXX_OBJECTS) $(C_OBJECTS) $(CC_OBJ)

#targets
all: GEM

GEM: $(CXX_OBJECTS) $(C_OBJECTS) $(CC_OBJ)
	$(CXX) $^ $(CXXFLAGS) $(BOOST_LIB_TD) $(BOOST_LIB) -o $@ $(DYN_LIBS)

%.o: %.c
	$(CC) -O2 -g -c $< -o $@

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -I$(BOOST_INC) -I$(Eigen_INC)

.PHONY: clean
clean:
	rm -f $(OBJS)
