CC = g++
#CFLAGS = -O3 -std=c++11 -Wall -fexceptions -fexpensive-optimizations -c
CFLAGS = -O0 -std=c++11 -Wall -fexceptions -c -g

all: particles psdistance

particles: main.o psystem.o intersections.o
	$(CC) main.o psystem.o intersections/intersections.o -o particles

psdistance: psdistance.o psystem.o intersections.o
	$(CC) psdistance.o psystem.o intersections/intersections.o -o psdistance

main.o: main.cpp psystem.h
	$(CC) $(CFLAGS) main.cpp -o main.o

psystem.o: psystem.cpp psystem.h intersections/intersections.h
	$(CC) $(CFLAGS) psystem.cpp -o psystem.o

intersections.o: intersections/intersections.cpp intersections/intersections.h
	$(CC) $(CFLAGS) intersections/intersections.cpp -o intersections/intersections.o

psdistance.o: psdistance.cpp psystem.h
	$(CC) $(CFLAGS) psdistance.cpp -o psdistance.o	

clean:
	rm -rf *.o particles psdistance intersections/intersections.o
