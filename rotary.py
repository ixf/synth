

import RPi.GPIO as GPIO
import time

pa = 18
pb = 23

value = 0
ok_states = { (0,0,0,1), (0,1,1,1), (1,1,1,0), (1,0,0,0) }

GPIO.setmode(GPIO.BCM)

def cback(channel):
	global last_a, last_b, value
	a, b = GPIO.input(pa), GPIO.input(pb)
	left = (last_a, last_b, a, b)
	right = (a, b, last_a, last_b)
	last_a, last_b = a, b
	if left in ok_states:
		value -= 1
		print(value)
	elif right in ok_states:
		value += 1
		print(value)
	else:
		print('?')
	

GPIO.setup(pa, GPIO.IN, pull_up_down = GPIO.PUD_UP)
GPIO.setup(pb, GPIO.IN, pull_up_down = GPIO.PUD_UP)

last_a, last_b = GPIO.input(pa), GPIO.input(pb)

GPIO.add_event_detect(pa, GPIO.BOTH, callback=cback)
GPIO.add_event_detect(pb, GPIO.BOTH, callback=cback)

time.sleep(60)

