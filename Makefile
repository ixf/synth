FLAGI=-lasound -lm -lncurses -I./PIGPIO/ -L./PIGPIO/ -l:libpigpio.so

all: e.c
	gcc e.c -o e $(FLAGI)

o: e.c
	gcc e.c -o e $(FLAGI) -O3

g: e.c
	gcc e.c -o e $(FLAGI) -g

rotary: rotary.c
	gcc rotary.c -o rotary $(FLAGI)
