CXX=g++

FLAGS= -lz

all: main

bfgs.o : bfgs.h bfgs.cpp
	$(CXX) -c bfgs.cpp $(FLAGS)

analysisfunctions.o: analysisfunctions.h analysisfunctions.cpp
	$(CXX) -c analysisfunctions.cpp $(FLAGS)

filehandlingfunctions.o: filehandlingfunctions.h filehandlingfunctions.cpp
	$(CXX) -c filehandlingfunctions.cpp $(FLAGS)

localngsrelate.o: localngsrelate.cpp
	$(CXX) -c localngsrelate.cpp $(FLAGS)

main: analysisfunctions.o localngsrelate.o filehandlingfunctions.o bfgs.o
	$(CXX) localngsrelate.o analysisfunctions.o filehandlingfunctions.o bfgs.o  $(FLAGS) -o localngsrelate

clean: 
	rm -f *.o localngsrelate
