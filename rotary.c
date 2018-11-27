

#include <stdlib.h>
#include <stdio.h>
#include <pigpio.h>

volatile int value = 0;

int last_gpio = 18;
int last_a = 0;
int last_b = 0;

char state = 0x00;

void callback(int GPIO, int level, unsigned int tick){
	printf("%d %d %d\n", value, GPIO, level);
	state <<= 2;
	if(level)
		if(GPIO == 18)
			state |= 2;
		else
			state |= 1;
	state &= 0xF;
	// 0001
	// 0111
	// 1110
	// 1000

	// 0100
	// 1101
	// 1011
	// 0010
	if(state == 1 || state == 7 || state == 14 || state == 8)
		value += 1;
	else if(state == 4 || state == 13 || state == 2 || state == 11)
		value -= 1;
}

int main(){

	if (gpioInitialise()<0){
		printf("init failed\n");
		return 1;
	}

	int pa = 18;
	int pb = 23;

	gpioSetMode(pa, PI_INPUT);
	gpioSetMode(pb, PI_INPUT);

	gpioSetPullUpDown(pa, PI_PUD_UP);
	gpioSetPullUpDown(pb, PI_PUD_UP);

	gpioSetAlertFunc(pa, callback);
	gpioSetAlertFunc(pb, callback);

	while(value < 200){}

	gpioSetAlertFunc(pa, 0);
	gpioSetAlertFunc(pb, 0);

	gpioTerminate();

}


