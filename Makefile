CXXFLAGS = -std=c++11 -Wall

.PHONY: all clean

all: plotfill2

plotfill2: main.o libclipper.a
	g++ -o plotfill2 main.o -lclipper

main.o: main.cpp
	g++ $(CXXFLAGS) -c -o main.o main.cpp

libclipper.a: clipper.cpp
	g++ -c -o clipper.o clipper.cpp && ar -rs libclipper.a clipper.o && rm clipper.o
