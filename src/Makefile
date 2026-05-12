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

reference: compressible.cpp compressible.hpp
	$(CXX) $(OPTFLAGS) -DNTEST=200 -o compressible_ref compressible.cpp
	./compressible_ref
	for f in rho u v et; do cp $${f}.csv $${f}_ref.csv; done

test: compressible.cpp compressible.hpp
	$(CXX) $(OPTFLAGS) -DNTEST=200 -o compressible_test compressible.cpp
	./compressible_test
	python3 -c "\
import numpy as np, sys; \
errs = {q: np.max(np.abs(np.loadtxt(q+'_ref.csv',delimiter=',') - np.loadtxt(q+'.csv',delimiter=','))) for q in ['rho','u','v','et']}; \
[print(f'{q}: {e:.2e}') for q,e in errs.items()]; \
sys.exit(0 if all(e < 1e-10 for e in errs.values()) else 1)"

.PHONY: clean reference test
clean:
	rm *.out
