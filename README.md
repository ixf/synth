

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

} Sound;

typedef struct {

  //czas uderzenia, czas wypuszczenia klawisza, klawisz np. A4

} Note;

//globale:

Sound currentSound;
Note 
Note currentNotes[];
Macro macros[];

combine_sounds(){

  //res = oblicz coś na podstawie aktualnego głównego dźwięku i jego trwających nut
        //i dodaj jakieś bajty z każdego z makr

}
```




