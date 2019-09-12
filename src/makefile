CC=g++
OBJECTS= MatrixUtils.o ReadParameters.o GEM.o
LIBS= -lboost_iostreams -lz -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -fopenmp -ldl -lpthread -lm

#targets
all: GEM
GEM: $(OBJECTS)
	$(CC) -o GEM $(OBJECTS) $(LIBS)
 
MatrixUtils.o: MatrixUtils.cpp
	$(CC) -c MatrixUtils.cpp

ReadParameters.o: ReadParameters.cpp
	$(CC) -c -std=c++11 ReadParameters.cpp

GEM.o: GEM.cpp
	$(CC) -std=c++11 -c -fopenmp GEM.cpp

#clean
clean:
	rm -f GEM $(OBJECTS)
