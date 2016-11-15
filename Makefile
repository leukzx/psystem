CC = g++
CFLAGS = -O3 -std=c++11 -Wall -fexceptions -fexpensive-optimizations -c

all: particles psdistance

particles: main.o psystem.o
	$(CC) main.o psystem.o -o particles

main.o: main.cpp psystem.h
	$(CC) $(CFLAGS) main.cpp -o main.o

psystem.o: psystem.cpp psystem.h intersections/intersections.h
	$(CC) $(CFLAGS) -c psystem.cpp -o psystem.o

intersections.o: intersections/intersections.cpp intersections/intersections.h
	$(CC) $(CFLAGS) -c intersections/intersections.cpp -o intersections/intersections.o

psdistance.o: psdistance.cpp psystem.h
	$(CC) $(CFLAGS) -c psdistance.cpp -o psdistance.o	

psdistance: psdistance.o psystem.o
	$(CC) psdistance.o psystem.o -o psdistance

clean:
	rm -rf *.o particles psdistance
