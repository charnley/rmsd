
FLAGS=-lblas -llapack -llapacke

all: matrix.o rmsd.o
	g++ -o calculate_rmsd ${FLAGS} *.o

calculate_rmsd: main.o kabsch.h
	g++ -o calculate_rmsd ${FLAGS} main.o

# g++ -o $@ ${FLAGS} $<
%.o: %.cpp
	g++ -c -Wall $<

python:
	c++ -O3 -Wall -shared -std=c++14 -fPIC `python3 -m pybind11 --includes` matrix.cpp binding.cpp -o example`python3-config --extension-suffix`

python_binding:
	c++ -O3 -Wall -shared -std=c++14 -fPIC `python3 -m pybind11 --includes` python_binding.cpp -o example`python3-config --extension-suffix`


clean:
	rm *.o *.so calculate_rmsd

