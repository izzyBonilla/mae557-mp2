# C++ Compiler Choice
CXX = g++

# Will look into adding parallelism later
OPTFLAGS = -O3
DEBUGFLAGS = -g

default: compressible.cpp compressible.hpp
	$(CXX) $(OPTFLAGS) -o compressible compressible.cpp

par: compressible.cpp compressible.hpp
	$(CXX) $(OPTFLAGS) -fopenmp -o compressible compressible.cpp

debug: compressible.cpp compressible.hpp
	$(CXX) $(DEBUGFLAGS) -o  compressible compressible.cpp

.PHONY: clean
clean:
	rm *.out
