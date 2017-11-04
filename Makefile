.PHONY: all clean

all: plotfill2

plotfill2: $(OBJS) clipper.a

clipper.a: clipper.cpp
	g++ -Wall -c -o clipper.o clipper.cpp && ar -rs clipper.a clipper.o && rm clipper.o
