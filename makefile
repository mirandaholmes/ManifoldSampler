CC=g++ -std=gnu++11 

#PROGRAM = examplepolymer
#PROGRAM = example1
PROGRAM = example2

OPTFLAGS =   -O3
EXECUTABLE=./$(PROGRAM) 

CPP_HEADERS = Equations.hpp Framework.hpp ManifoldSampler.hpp

all: $(EXECUTABLE) run 

$(PROGRAM): $(PROGRAM).o ManifoldSampler.o Framework.o
	$(CC) $(OPTFLAGS) -o  $(PROGRAM) $^

$(PROGRAM).o: $(PROGRAM).cpp $(CPP_HEADERS) makefile
	$(CC) -c $(OPTFLAGS) $(PROGRAM).cpp 

ManifoldSampler.o:  ManifoldSampler.cpp ManifoldSampler.hpp Equations.hpp 
	$(CC) -c $(OPTFLAGS) ManifoldSampler.cpp

Framework.o:  Framework.cpp Framework.hpp Equations.hpp
	$(CC) -c $(OPTFLAGS) Framework.cpp

run: $(EXECUTABLE)
	$(EXECUTABLE) $(ARGS) 

clean:
	rm -f *.o  