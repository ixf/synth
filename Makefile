FLAGI=-lasound -lm -lncurses

all: alsa_stuff.c e.c
	gcc $^ -o e $(FLAGI)

pi: alsa_stuff.c e.c
	gcc $^ -o e $(FLAGI) -I./PIGPIO/ -L./PIGPIO/ -l:libpigpio.so -DMAKE_PI -O3

o: alsa_stuff.c e.c
	gcc $^ -o e $(FLAGI) -O3
	
g: alsa_stuff.c e.c
	gcc $^ -o e $(FLAGI) -g

