
#include <stdlib.h>
#include <stdio.h>

#ifdef MAKE_PI
#include <pigpio.h>
#endif

#include <sys/time.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/prctl.h>

#include <string.h>
#include <signal.h>
#include <sched.h>
#include <errno.h>
#include <unistd.h>
#include <getopt.h>
#include <alsa/asoundlib.h>
#include <math.h>

#include <ncurses.h>

#include <fcntl.h>
#include <dirent.h>
#include <linux/input.h>

#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/select.h>
#include <termios.h>

long long get_time() {
    struct timeval te; 
    gettimeofday(&te, NULL);
    long long microseconds = te.tv_sec*1000000LL + te.tv_usec;
    return microseconds;
}
 
#include "alsa_stuff.h"
#include "main.h"

FILE* logfile;

// indeks elementu z shared_values który ma ma modyfikować n-ty enkoder
// ustawiana przez główny proces, używa fork
static int* shared_indexes;

typedef struct {
  double min;
  double max;
  double step;
} range;

range value_ranges[6] = {
  {0, 8000, 100},
  {0, 8000, 100},
  {0, 5, 0.01},
  {0, 5, 0.01},
  {0, 1, 0.01},
  {0, 5, 0.01}
};

// tablica zmiennych ( przyjęliśmy że z zakresu 0-255 )
// wartości zmienia fork, odczytuje główny
static double* shared_values;

#ifdef MAKE_PI
typedef struct {
  int last_a;
  int last_b;
  char state;
} rotary_state;

rotary_state rotary_states[6];
int gpio_to_n[40]; // 0 do 6
int gpio_to_ab[40]; // 0 lub 1

void rot_callback(int GPIO, int level, unsigned int tick){

  // n-ty enkoder, wejście a lub b
  int n = gpio_to_n[GPIO];
  int ab = gpio_to_ab[GPIO];

  rotary_state* rs = &(rotary_states[n]);

  printf("rot callback!\n");

  rs->state <<= 2;
  if(ab == 0){
    if(level)
      rs->state |= 0b0010;
    if(rs->last_b)
      rs->state |= 0b0001;
  } else {
    if(rs->last_a)
      rs->state |= 0b0010;
    if(level)
      rs->state |= 0b0001;
  }
  rs->state &= 0xF;

  // 0001
  // 0111
  // 1110
  // 1000

  // 0100
  // 1101
  // 1011
  // 0010
  
  int shared_index = shared_indexes[n];
  range r = value_ranges[shared_index];

  printf("< sv[%d] := %lf\r\n", shared_index, shared_values[shared_index]);
  if(rs->state == 1 || rs->state == 7 || rs->state == 14 || rs->state == 8){
    double now = shared_values[shared_index] += r.step;
    printf("> sv[%d] := %lf\r\n", shared_index, now);
  } else if(rs->state == 4 || rs->state == 13 || rs->state == 2 || rs->state == 11){
    double now = shared_values[shared_index] -= r.step;
    printf("sv[%d] := %lf\n", shared_index, now);
  }

  if(shared_values[shared_index] > r.max)
    shared_values[shared_index] = r.max;
  else if (shared_values[shared_index] < r.min)
    shared_values[shared_index] = r.min;
}

void setup_rotary_encoder(int pa, int pb, int n){
    gpio_to_ab[pa] = 0;
    gpio_to_ab[pb] = 1;
    gpio_to_n[pa] = n;
    gpio_to_n[pb] = n;

    gpioSetMode(pa, PI_INPUT);
    gpioSetMode(pb, PI_INPUT);

    gpioSetPullUpDown(pa, PI_PUD_UP);
    gpioSetPullUpDown(pb, PI_PUD_UP);

    gpioSetAlertFunc(pa, rot_callback);
    gpioSetAlertFunc(pb, rot_callback);

    /* gpioSetAlertFunc(pa, 0); */
    /* gpioSetAlertFunc(pb, 0); */
    /* gpioTerminate(); */
}

void child_setup(){
  setup_rotary_encoder(18,23,0);
  setup_rotary_encoder(14,15,1);

  shared_indexes[0] = 0;
  shared_indexes[1] = 2;
}
#endif 

// **************************************** FALE

double sin_wave(double x){
  return sinf(x);
}

double step3_wave(double x){
  return ceil(3*x/(2*M_PI))-2;
}

double step8_wave(double x){
  return ceil(8*x/(2*M_PI) - 4.5)/3.5;
}

double triangle_wave(double x){
  if( x < 0.5 * M_PI )
    return x*2/M_PI;
  else if( x < 1.5 * M_PI )
    return (M_PI-x)*2/M_PI;
  else
    return x*2/M_PI-4.0;
}

double sawup_wave(double x){
  return 1.0 - x/M_PI;
}

double sawdown_wave(double x){
  return x/M_PI - 1.0;
}

double square_wave(double x){
  if(x > M_PI)
    return 1.0;
  else
    return -1.0;
}

typedef double (*wave_ptr)(double x);
typedef struct {
  char* name;
  wave_ptr fun;
} wave;

wave waves[7] = {
  {"sinus", sin_wave },
  {"sawup", sawup_wave },
  {"sawdown", sawdown_wave },
  {"schodki 3", step3_wave} ,
  {"schodki 8", step8_wave },
  {"trójkątna", triangle_wave },
  {"prostokątna", square_wave },
};

double weight[2] = { 0.5, 0.5 };
int wave1_index = 0;
int wave2_index = 1;

//**************************************** FILTR

typedef struct {
  double* flow;
  double* fhigh;

  double last_flow;
  double last_fhigh;

  double phi;
  double lambda;
  int size;

  double* samples;
  double* params;
} BandpassFilter;

void refresh_bpf(BandpassFilter* filter){
  if( filter->last_flow != *(filter->flow) ||
   filter->last_fhigh != *(filter->fhigh) ){
    // refresh
    
    filter->last_flow = *(filter->flow);
    filter->last_fhigh = *(filter->fhigh);


    filter->lambda = M_PI * filter->last_flow / (44100/2);
    filter->phi = M_PI * filter->last_fhigh / (44100/2);

    if(filter->samples != NULL){
      free(filter->samples);
      free(filter->params);
    }

    filter->samples = malloc(sizeof(double) * filter->size);
    filter->params = malloc(sizeof(double) * filter->size);

    for(int i = 0; i < filter->size; i++){
      filter->samples[i] = 0;

      double dist = i - (filter->size - 1.0) / 2.0;
      if( dist == 0.0 ){
	filter->params[i] = (filter->phi - filter->lambda) / M_PI;
      } else {
	filter->params[i] = ( sin( dist * filter->phi ) - sin( dist * filter->lambda )) / (dist * M_PI);
      }
    }
  }
}

void init_bpf(BandpassFilter* filter, double* flow, double* fhigh, int size){
  filter->flow = flow;
  filter->fhigh = fhigh;
  filter->size = size;

  refresh_bpf(filter);
}

double get_filtered_sample(BandpassFilter* f, double in){

  for(int i = f->size - 1; i > 0; i--){
    f->samples[i] = f->samples[i-1];
  }

  f->samples[0] = in;

  double result = 0;
  for(int i = 0; i < f->size; i++){
    result += f->samples[i] * f->params[i];
  }

  return result;
}


typedef struct {
	long long attack;
	long long release;
	double freq;
	double step;
	double phase;

	bool active;
} Note;

//**************************************** ADSR

double bad_adsr(Note *n, long long now){
	if(n->release == -1){
		return 1.0;
	}
	double seconds = 0.2;
	double r = fmin(1.0, fmax(0.0, 1.0 - ((double)(now - n->release))/CLOCKS_PER_SEC/(seconds)));
	if ( r == 0.0 )
		n->active = false;
	return r;
}

// parametry:
// A -- ATTACK -- maks czas rośnięcia od 0.0 do 1.0 kiedy klawisz jest wciśnięty
// D -- RELEASE -- czas wyciszania kiedy klawisz jest nadal wciśnięty
// S -- SUSTAIN -- głośność w jakiej nuta się utrzyma kiedy klawisz będzie przytrzymany
// R -- RELEASE -- czas wyciszania po puszczeniu klawisza

typedef struct {
  double* a;
  double* d;
  double* s;
  double* r;
} adsr_params;

adsr_params main_adsr_params;

double lin_adsr(Note *n, long long now){
  double x = (now - n->attack)/((double)CLOCKS_PER_SEC);
  double a = *(main_adsr_params.a), d = *(main_adsr_params.d),
	 s = *(main_adsr_params.s), r = *(main_adsr_params.s);

  //fprintf(logfile, "%lf %lf %lf %lf\n", a,d,s,r);
  if( n->release == -1 ){

    if( x > a + d ){
      // etap sustain
      return s;
    } else if( x < a ) {
      // etap attack
      return x / a;
    } else {
      // delay
      return 1-(1.0-s)*(x-a)/d;
    }
    
  } else {
    double x0 = (n->release - n->attack)/((double)CLOCKS_PER_SEC);
    double x2 = (now - n->release)/((double)CLOCKS_PER_SEC);

    if( x2 > r ){
      n->active = false;
      return 0.0;
    } else if( x0 >= a+d ){
      // release po sustain
      
      return s * (r-x2)/r;
    } else {
      // release przed sustain
      //
      // aby było gładko itd należy przez czas r wygłuszać dźwięk od
      // poprzedniej wartości do 0
      // poprzednią wartość obliczamy jakby n->release było równe -1
      // a naszym nowym x jest czas od released do teraz
      double y;

      if( x0 < a ) {
	y = x0 / a;
      } else {
	y = 1-(1.0-s)*(x0-a)/d;
      }


      return y * (r-x2)/r;

    }
  }
}

double freq_calc(int n){ // A4 = 49 -> 440Hz
	return 440.0 * pow(2.0, (n-49)/12.0);
}

char piano_keys[255] = { 0 };
void init_piano_keys(int starting, int* keys, int count){
	// przyklad: 40, "awsed"
	// piano_keys['a'] = 40
	// piano_keys['w'] = 41
	// ...

	for(int i = 0; i < count; i++){
		piano_keys[keys[i]] = starting+i;
	}
}

double (*main_adsr)(Note* note, long long now);
BandpassFilter main_filter;
Note all_notes[104];

static double max_phase = 2. * M_PI;

void control_loop(){

  for(int i = 0; i < 104; i++){
    all_notes[i].active = false;
    double freq = freq_calc(i);
    all_notes[i].freq = freq;
    all_notes[i].step = max_phase*(freq)/44100.0;
  }


  struct input_event ev[1];
  int fd, rd, value, code, size = sizeof (struct input_event);
  char name[256] = "Unknown";
  char *device = NULL;

  if ((getuid ()) != 0)
    printf ("You are not root! This may not work...n/");

#if MAKE_PI
  device = "/dev/input/event0";
#else
  device = "/dev/input/event4";
#endif
  //Open Device
  if ((fd = open (device, O_RDONLY)) == -1)
    printf ("%s is not a vaild device.n", device);
  ioctl (fd, EVIOCGNAME (sizeof (name)), name);
  printf ("IO: Reading From : %s (%s)\n", device, name);
  fflush(stdout); 


  while (1){
    if ((rd = read (fd, ev, size )) < size){
      perror("Error reading");  
      abort();
    }
    value = ev[0].value;
    code = ev[0].code;
    if (value !=2 && ev[0].type == 1){
      //got char
      /* printf ("Code[%d] %d \r\n", code, value); */
      /* fflush(stdout); */
      /* printf ("Code = %d \r\n",code); */

      if( value == 1 && code == 59 ){
	printf("F1 set 0: \r\n");
	//todo
	continue;
      } else if (value == 1 && code == 60){
	printf("F2 set 1: \r\n");
	continue;
      } else if (value == 1 && code == 61){
	wave1_index = (wave1_index+1)%7;
	printf("F3 fala 1 to teraz: %s\r\n", waves[wave1_index].name);
	continue;
      } else if (value == 1 && code == 62){
	wave2_index = (wave2_index+1)%7;
	printf("F4 fala 2 to teraz: %s\r\n", waves[wave2_index].name);
	continue;
      }

      int note = piano_keys[code];
      if(value == 1){
	//PRESS
	all_notes[note].phase = 0.0;
	all_notes[note].active = true;
	all_notes[note].attack = get_time();
	all_notes[note].release = -1;
      }else{
	//RELEASE
	all_notes[note].release = get_time();
      }
    }
  }  
  exit(EXIT_SUCCESS);
}

void before_samples_loop(){

  refresh_bpf(&main_filter);

}

double get_new_sample(){
  int format_bits = snd_pcm_format_width(format);
  unsigned int maxval = (1 << (format_bits - 1)) - 1;

      int res = 0;

      long long clock_now = get_time();
      for(int i = 28; i <= 68; i++){
	Note* n = &(all_notes[i]);

	if(! n->active )
	  continue;

	double adsr_val = main_adsr(n, clock_now);

	// suma z dwóch fal:
	wave w1 = waves[wave1_index];
	wave w2 = waves[wave2_index];
	double osc_total = w1.fun(n->phase) * weight[0] + w2.fun(n->phase) * weight[1];

	// efekt ADSR
	osc_total = osc_total * adsr_val * maxval;

	// todo maxval przenieść?

	n->phase += n->step;
	if (n->phase >= max_phase){
	  n->phase -= max_phase;
	}

	res += osc_total;
      }

      res /= 10.0;

      // teraz res jest wartością z sumy wszystkich Note
      // wrzucamy do filtra bandpass:

      res = get_filtered_sample(&main_filter, res);

      return res;
}



int main(int argc, char *argv[]) {

  logfile = fopen("logfile", "w");

  shared_indexes = mmap(NULL, 6*sizeof(int), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
  shared_values = mmap(NULL, 32*sizeof(double), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

  shared_values[0] = 3000; // dolna granica filtra
  shared_values[1] = 4000; // górna granica filtra

  shared_values[2] = 0.05; // a
  shared_values[3] = 0.05; // d
  shared_values[4] = 0.1; // s
  shared_values[5] = 0.0; // r

  shared_values [6] = 0.5; // master volume
  shared_values [7] = 0.5; // waga fali 1
  shared_values [8] = 0.5; // waga fali 2

  main_adsr_params.a = &(shared_values[2]);
  main_adsr_params.d = &(shared_values[3]);
  main_adsr_params.s = &(shared_values[4]);
  main_adsr_params.r = &(shared_values[5]);

  main_adsr = lin_adsr;
  init_bpf( &(main_filter), &(shared_values[0]), &(shared_values[1]), 17);

#ifdef MAKE_PI
  int child;

  if((child = fork()) == 0){
    prctl(PR_SET_PDEATHSIG, SIGHUP);
    if (gpioInitialise()<0){
      printf("init failed\n");
      return 1;
    }

    child_setup();
    printf("while1 started\r\n");
    while(1){}

    exit(EXIT_SUCCESS);
  }
#endif


  initscr();

  int top_row[] = { 16, 3,  17,  4, 18, 19,  6, 20,  7, 21,  8, 22 };
  int bot_row[] = { 44, 31, 45, 32, 46, 47, 34, 48, 35, 49, 36, 50 };
  init_piano_keys(52, top_row, sizeof(top_row)/sizeof(top_row[0]));
  init_piano_keys(40, bot_row, sizeof(bot_row)/sizeof(bot_row[0]));

  int alsa_status = init_alsa(argc, argv);
  if( alsa_status != 0 ){
    printf("init_alsa zwróciło %d\n", alsa_status);
    return alsa_status;
  }

  int err = alsa_loop();

  if (err < 0)
    printf("Transfer failed: %s\n", snd_strerror(err));

  alsa_close();

#ifdef MAKE_PI
  kill(child, 9);
#endif
  wait(NULL);
  munmap(shared_indexes, sizeof(int)*6);
  munmap(shared_values, sizeof(double)*32);

  endwin();
  return 0;
}
