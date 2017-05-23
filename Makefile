CPP=g++
CPPFLAGS=-O3 -std=c++14 -fopenmp -lpthread
#CPP=icpc
#CPPFLAGS=-O3 -qopenmp -std=c++14

all:
	${CPP} ${CPPFLAGS} example2.cpp -o example2
clean:
	rm -f Makefile~ example2.cpp~
	rm -f example2

prof:
	valgrind --dump-instr=yes --tool=callgrind ./example2 1 10

#
