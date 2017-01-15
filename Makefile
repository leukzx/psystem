CC = g++
CFLAGS = -c -O3 -fopenmp -std=c++11 -Wall -fexceptions -fexpensive-optimizations
BUILDDIR = ./build
LIB_DIR = ./ 
INCLUDE_DIR = ./build
LIBS = -lconfig++ -lgomp -lOpenCL
 
all: psystem psdistance

psystem: main.o psystem.o intersections.o
	$(CC) -I$(INCLUDE_DIR) -L$(LIB_DIR) -o $(BUILDDIR)/psystem $(BUILDDIR)/main.o $(BUILDDIR)/psystem.o $(BUILDDIR)/intersections.o $(LIBS)

psdistance: psdistance.o psystem.o intersections.o
	$(CC) -I$(INCLUDE_DIR) -L$(LIB_DIR) -o $(BUILDDIR)/psdistance $(BUILDDIR)/psdistance.o $(BUILDDIR)/psystem.o $(BUILDDIR)/intersections.o $(LIBS)

main.o: main.cpp psystem.h
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) main.cpp -o $(BUILDDIR)/main.o

psystem.o: psystem.cpp psystem.h intersections/intersections.h
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) psystem.cpp -o $(BUILDDIR)/psystem.o

intersections.o: intersections/intersections.cpp intersections/intersections.h
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) intersections/intersections.cpp -o $(BUILDDIR)/intersections.o

psdistance.o: psdistance.cpp psystem.h
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) psdistance.cpp -o $(BUILDDIR)/psdistance.o	

clean:
	rm -rf  $(BUILDDIR)/*.o $(BUILDDIR)/psdistance $(BUILDDIR)/psystem
