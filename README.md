
## Syntezator na raspberry pi

Wykorzystujemy: ALSA i PIGPIO.

Etapy powstawania dźwięku:
 * na raz gra wiele nut o różnych częstotliwościach
 * mamy dwie funkcje bazowe ( wybrane z 7 ) o dwóch wagach w oparciu o fazę poszczególnej nuty
 * mnożymy powyższe razy głośność nuty ( ADSR )
 * sumujemy wartości wszystkich nut
 * wynik wrzucamy do filtra bandpass

Parametry ADSR i filtra można na bierząco zmieniać poprzez pokrętła
