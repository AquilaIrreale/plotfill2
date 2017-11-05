CXXFLAGS = -std=c++11 -Wall
LDFLAGS  = -L. -lclipper -lpng

.PHONY: all clean

all: plotfill2

plotfill2: main.o libclipper.a
	g++ $(CXXFLAGS) -o plotfill2 main.o $(LDFLAGS)

main.o: main.cpp
	g++ $(CXXFLAGS) -c -o main.o main.cpp

libclipper.a: clipper.cpp
	g++ -c -o clipper.o clipper.cpp && ar -rs libclipper.a clipper.o

clean:
	@rm -rvf plotfill2 *.a *.o
