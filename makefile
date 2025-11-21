CXX=clang++
OPT_LEVEL=3
STANDART=c++98

OMP_CXXFLAGS=-O$(OPT_LEVEL) -march=native -mavx2 -fopenmp -ggdb -Wall -Wextra -pedantic -std=$(STANDART)
CXXFLAGS=-O$(OPT_LEVEL) -march=native -mavx2 -ggdb -Wall -Wextra -pedantic -std=$(STANDART)

all: visual main main_omp

pgo: main_pgo main_omp_pgo

visual: STANDART=c++17
visual: OPT_LEVEL=0
visual: visual.cpp main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< -lraylib

main: main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

main_omp: main.cpp
	$(CXX) $(OMP_CXXFLAGS) -o $@ $<

main_pgo: main.cpp main.instr
	./main.instr
	llvm-profdata merge default.profraw -o main_pgo.profdata
	$(CXX) -fprofile-instr-use=main_pgo.profdata $(CXXFLAGS) $< -o $@

main.instr: OPT_LEVEL = 2
main.instr: main.cpp
	$(CXX) -fprofile-instr-generate $(CXXFLAGS) $< -o main.instr

main_omp_pgo: main.cpp main_omp.instr
	./main_omp.instr
	llvm-profdata merge default.profraw -o main_omp_pgo.profdata
	$(CXX) -fprofile-instr-use=main_omp_pgo.profdata $(OMP_CXXFLAGS) $< -o $@

main_omp.instr: OPT_LEVEL = 2
main_omp.instr: main.cpp
	$(CXX) -fprofile-instr-generate $(OMP_CXXFLAGS) $< -o main_omp.instr
