CXX=clang++
OPT_LEVEL=3
OMP_CXXFLAGS=-O$(OPT_LEVEL) -march=native -mavx2 -fopenmp -ggdb -Wall -Wextra -pedantic -std=c++98
CXXFLAGS=-O$(OPT_LEVEL) -march=native -mavx2 -ggdb -Wall -Wextra -pedantic -std=c++98

all: main main_omp

main: main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

main_omp: main.cpp
	$(CXX) $(OMP_CXXFLAGS) -o $@ $<
