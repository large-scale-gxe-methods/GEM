CC=g++
OBJECTS= MatrixUtils.o ReadParameters.o TimeUtils.o ReadBGEN.o GEM.o
LIBS= -I../thirdparty/libdivide -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -ldl -lpthread -lm -lboost_program_options -lboost_thread -lboost_system -lboost_filesystem -lzstd -L../thirdparty/libdeflate -ldeflate -Wl,-rpath=$(PWD)/../thirdparty/libdeflate



#targets
all: GEM


GEM: $(OBJECTS)
	$(MAKE) -C ../thirdparty/libdeflate
	$(CC) -o GEM $(OBJECTS) $(LIBS)

MatrixUtils.o: MatrixUtils.cpp
	$(CC) -c -O2 MatrixUtils.cpp

ReadParameters.o: ReadParameters.cpp
	$(CC) -c -std=c++11 ReadParameters.cpp

TimeUtils.o: TimeUtils.cpp
	$(CC) -c -std=c++11 TimeUtils.cpp

ReadBGEN.o: ReadBGEN.cpp
	$(CC) -c -std=c++11 -O2 -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0 ReadBGEN.cpp

GEM.o: GEM.cpp
	$(CC) -std=c++11 -c -O2 GEM.cpp

#clean
clean:
	rm -f GEM $(OBJECTS)
