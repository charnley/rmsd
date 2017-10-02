
FLAGS=-lblas -llapack -llapacke

all: rmsd.o

%.o: %.cpp
	g++ -o $@ ${FLAGS} $<

