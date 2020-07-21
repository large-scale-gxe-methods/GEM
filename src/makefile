CC=gcc
CXX=g++
CXXFLAGS = -std=c++11 -O2 -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0
LIBS= -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -ldl -lpthread -lm -lboost_program_options -lboost_thread -lboost_system -lboost_filesystem

SOURCES= $(wildcard ../thirdparty/zstd-1.4.5/lib/common/*.c) $(wildcard ../thirdparty/zstd-1.4.5/lib/compress/*.c) $(wildcard ../thirdparty/zstd-1.4.5/lib/decompress/*.c) $(wildcard ../thirdparty/zstd-1.4.5/lib/dictBuilder/*.c) $(wildcard ../thirdparty/libdeflate-1.6/lib/*.c) $(wildcard ../thirdparty/libdeflate-1.6/lib/x86/*.c) $(wildcard ../thirdparty/libdeflate-1.6/lib/arm/*.c)
OBJS= MatrixUtils.o ReadParameters.o TimeUtils.o ReadBGEN.o GEM.o $(SOURCES:.c=.o) 

#targets
all: GEM


GEM: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIBS)


#clean
	rm -f $(OBJS)

