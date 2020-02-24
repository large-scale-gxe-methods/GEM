CC=g++
OBJECTS= MatrixUtils.o ReadParameters.o TimeUtils.o ReadBGEN.o GEM.o
LIBS= -lboost_iostreams -lz -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -ldl -fopenmp -lm -lboost_program_options -lboost_thread -lboost_system -lboost_filesystem

OMP_NUM_THREADS=1
export OMP_NUM_THREADS

#targets
all: GEM
GEM: $(OBJECTS)
	$(CC) -o GEM $(OBJECTS) $(LIBS)
 
MatrixUtils.o: MatrixUtils.cpp
	$(CC) -c MatrixUtils.cpp

ReadParameters.o: ReadParameters.cpp
	$(CC) -c -std=c++11 ReadParameters.cpp

TimeUtils.o: TimeUtils.cpp
	$(CC) -c -std=c++11 TimeUtils.cpp

ReadBGEN.o: ReadBGEN.cpp
	$(CC) -c -std=c++11 ReadBGEN.cpp

GEM.o: GEM.cpp
	$(CC) -std=c++11 -c -fopenmp GEM.cpp

#clean
clean:
	rm -f GEM $(OBJECTS)
