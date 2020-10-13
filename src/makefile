CC=gcc -std=gnu99 # Do not change this for libdeflate package.
CXX=g++
CXXFLAGS = -std=c++11 -O2 -g -Wall

LIBS= -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -ldl -lpthread -lm -lboost_program_options -lboost_thread -lboost_system -lboost_filesystem

C_SOURCES= $(wildcard ../thirdparty/zstd-1.4.5/lib/common/*.c) $(wildcard ../thirdparty/zstd-1.4.5/lib/compress/*.c) $(wildcard ../thirdparty/zstd-1.4.5/lib/decompress/*.c) $(wildcard ../thirdparty/zstd-1.4.5/lib/dictBuilder/*.c) $(wildcard ../thirdparty/libdeflate-1.6/lib/*.c) $(wildcard ../thirdparty/libdeflate-1.6/lib/x86/*.c) $(wildcard ../thirdparty/libdeflate-1.6/lib/arm/*.c)
CXX_SOURCES= MatrixUtils.cpp ReadParameters.cpp TimeUtils.cpp ReadBGEN.cpp ReadPGEN.cpp ReadBed.cpp GEM.cpp $(wildcard ../thirdparty/plink-2.0/*.cpp)
CXX_OBJECTS= $(CXX_SOURCES:.cpp=.o)
C_OBJECTS= $(C_SOURCES:.c=.o)
OBJS= $(CXX_OBJECTS) $(C_OBJECTS)

#targets
all: GEM

GEM: $(CXX_OBJECTS) $(C_OBJECTS)
	$(CXX) $^ $(CXXFLAGS) -o $@ $(LIBS)

%.o: %.c
	$(CC) -O2 -g -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJS)
