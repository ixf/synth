

dla każdego używanego dźwięku trzymamy strukturę

struktura zawiera parametry dźwięku typu rodzaj fali,
ADSR, harmoniczne, jakieś cuda

dynamiczna kompilacja/linkowanie wpisanego kodu w c jako biblioteka dynamiczna/dzielona dla beki

chyba że nasz dźwięk jest dźwiękiem właśnie granym z klawiatury


```c
typedef struct {

  Note notes[];
  byte* array;

  Sound s;

  time total_length;
  time phase;

} Makro;


typedef struct {

  //parametry z enkoderów
  // m.in.
  ADSR dla głośności
  filtry częstotliwości
  ADSR dla filtrów

} Sound;

typedef struct {

  //czas uderzenia, czas wypuszczenia klawisza, klawisz np. A4

} Note;

//globale:

Sound currentSound;
Note currentNotes[];
Macro macros[];

combine_sounds(){

  //res = oblicz coś na podstawie aktualnego głównego dźwięku i jego trwających nut
        //i dodaj jakieś bajty z każdego z makr

}
```

Kroki
jest funkcja bazowa
w jakiś stały prosty sposób dodajemy harmoniczne
usuwamy je filtrem
zmieniamy głośność w czasie z ADSR

TERAZ:
  - ogarnąć metody ( async direct itp ) i wyrzucić śmieci z kodu
  - terminalowy interfejs
  - kompilowanie/linkowanie funkcji bazowych
  - ADSR jak helm z 4 parametrami
    - funckja która do pewnego momentu rośnie a potem spada
    - punkt S podczas spadku
    - A liniowe, D i R śmieszne
  - filtry częstotliwości
  - sklejenie wszystkiego w combine_sounds

PÓŹNIEJ:
  - makra
  - LFO
  - efekty




