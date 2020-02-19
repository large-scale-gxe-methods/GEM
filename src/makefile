CC=g++
OBJECTS= MatrixUtils.o ReadParameters.o ReadBGEN.o TimeUtils.o GEM.o
LIBS= -lboost_iostreams -lz -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -ldl -lpthread -fopenmp -lm -lboost_program_options -lboost_thread -lboost_system -lboost_filesystem

#targets
all: GEM
GEM: $(OBJECTS)
	$(CC) -o GEM $(OBJECTS) $(LIBS)
 
MatrixUtils.o: MatrixUtils.cpp
	$(CC) -c MatrixUtils.cpp

ReadParameters.o: ReadParameters.cpp
	$(CC) -c -std=c++11 ReadParameters.cpp

ReadBGEN.o: ReadBGEN.cpp
	$(CC) -c -std=c++11 ReadBGEN.cpp

TimeUtils.o: TimeUtils.cpp
	$(CC) -C -std=c++11 TimeUtils.cpp

GEM.o: GEM.cpp
	$(CC) -std=c++11 -c -fopenmp GEM.cpp

#clean
clean:
	rm -f GEM $(OBJECTS)
