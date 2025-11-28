CC = gcc
CFLAGS = -Wall -O3
LIBS = -lm

all: main

main: main.o voronoi.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

main.o: main.c
	$(CC) $(CFLAGS) -o $@ -c $<

voronoi.o: voronoi.c
	$(CC) $(CFLAGS) -o $@ -c $<

test_voronoi.o: test_voronoi.c
	$(CC) $(CFLAGS) -o $@ -c $<

test_voronoi: test_voronoi.o voronoi.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

clean:
	rm *.o main
