#
# Makefile for ManifoldSampler. 
# Make sure you have Eigen installed in the compiler library path. 
#
#
#
#
# MHC Notes: 
#  on cims computers, run
#
#             module load gcc-9.2
#
#  first, to load the correct compiler
#
#

CC=g++ -std=gnu++11 

PROGRAM = examples_paper
#PROGRAM = run_poly
#PROGRAM = run_lattice
#PROGRAM = run_matrix
#PROGRAM = run_ngon
#PROGRAM = run_timing

#OPTFLAGS =   -O3  -I/home/holmes/libs   # on CIMS
OPTFLAGS =   -O3      # on PRO
EXECUTABLE=./$(PROGRAM) 

CPP_HEADERS = Equations.hpp ManifoldSampler.hpp Polygon.hpp Polymer.hpp Lattice.hpp Matrix1.hpp

all: $(EXECUTABLE) #run 

$(PROGRAM): $(PROGRAM).o ManifoldSampler.o 
	$(CC) $(OPTFLAGS) -o  $(PROGRAM) $^

$(PROGRAM).o: $(PROGRAM).cpp $(CPP_HEADERS) makefile
	$(CC) -c $(OPTFLAGS) $(PROGRAM).cpp 

ManifoldSampler.o:  ManifoldSampler.cpp ManifoldSampler.hpp Equations.hpp Polygon.hpp 
	$(CC) -c $(OPTFLAGS) ManifoldSampler.cpp

Framework.o:  Framework.cpp Equations.hpp Polygon.hpp
	$(CC) -c $(OPTFLAGS) Framework.cpp

run: $(EXECUTABLE)
	$(EXECUTABLE) $(ARGS) 

clean:
	rm -f *.o  