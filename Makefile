FLAGI=-lasound -lm -lncurses 

all: e.c
	gcc e.c -o e $(FLAGI)

o: e.c
	gcc e.c -o e $(FLAGI) -O3 --fast-math

g: e.c
	gcc e.c -o e $(FLAGI) -g

run: all
	./e -D plughw:0,0 -m async_direct

run01: all
	./e -D plughw:01,00 -m async_direct

