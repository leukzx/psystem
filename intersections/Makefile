CC = g++
CFLAGS = -O3 -std=c++11 -Wall -fexceptions -fexpensive-optimizations -c

all: intersections

intersections: main.o intersections.o
	$(CC) main.o intersections.o -o intersections

main.o: main.cpp intersections.h
	$(CC) $(CFLAGS) main.cpp -o main.o

intersections.o: intersections.cpp intersections.h
	$(CC) $(CFLAGS) -c intersections.cpp -o intersections.o

clean:
	rm -rf *.o intersections
