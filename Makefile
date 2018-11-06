FLAGI=-lasound -lm -lncurses 

all: e.c
	gcc e.c -o e $(FLAGI)

o: e.c
	gcc e.c -o e $(FLAGI) -O3

g: e.c
	gcc e.c -o e $(FLAGI) -g
