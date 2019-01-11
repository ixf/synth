FLAGI=-lasound -lm -lncurses

all: e.c
	gcc e.c -o e $(FLAGI)

pi: e.c
	gcc e.c -o e $(FLAGI) -I./PIGPIO/ -L./PIGPIO/ -l:libpigpio.so -DMAKE_PI

o: e.c
	gcc e.c -o e $(FLAGI) -O3

g: e.c
	gcc e.c -o e $(FLAGI) -g

rotary: rotary.c
	gcc rotary.c -o rotary $(FLAGI)
