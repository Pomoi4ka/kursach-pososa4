CXX=clang++
OPT_LEVEL=3
OMP_CXXFLAGS=-O$(OPT_LEVEL) -march=native -mavx2 -fopenmp -ggdb -Wall -Wextra -pedantic -std=c++98
CXXFLAGS=-O$(OPT_LEVEL) -march=native -mavx2 -ggdb -Wall -Wextra -pedantic -std=c++98

all: main main_omp

pgo: main_pgo main_omp_pgo

main: main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

main_omp: main.cpp
	$(CXX) $(OMP_CXXFLAGS) -o $@ $<

main_pgo: main.cpp main.instr
	./main.instr
	llvm-profdata merge default.profraw -o main_pgo.profdata
	$(CXX) -fprofile-instr-use=main_pgo.profdata $(CXXFLAGS) $< -o $@

main.instr: main.cpp
	$(CXX) -fprofile-instr-generate $(CXXFLAGS) $< -o main.instr

main_omp_pgo: main.cpp main_omp.instr
	./main_omp.instr
	llvm-profdata merge default.profraw -o main_omp_pgo.profdata
	$(CXX) -fprofile-instr-use=main_omp_pgo.profdata $(OMP_CXXFLAGS) $< -o $@

main_omp.instr: main.cpp
	$(CXX) -fprofile-instr-generate $(OMP_CXXFLAGS) $< -o main_omp.instr
