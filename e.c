
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
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/select.h>
#include <sys/time.h>
#include <termios.h>
 

// dzielone
static int *rot_state;

// dla forka:

int last_gpio = 18;
int last_a = 0;
int last_b = 0;

char state = 0x00;


// 0001
// 0111
// 1110
// 1000

// 0100
// 1101
// 1011
// 0010

void rot_callback(int GPIO, int level, unsigned int tick){
	//printf("%d %d %d\n", *rot_state, GPIO, level);
	state <<= 2;
	if(level)
		if(GPIO == 18)
			state |= 2;
		else
			state |= 1;
	state &= 0xF;

	if(state == 1 || state == 7 || state == 14 || state == 8)
		*rot_state += 1;
	else if(state == 4 || state == 13 || state == 2 || state == 11)
		*rot_state -= 1;
}

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
double weight[2] = { 0.5, 0.5 };
wave_ptr osc_wave1 = sin_wave;
wave_ptr osc_wave2 = sawup_wave;

typedef struct {
  double phi;
  double lambda;
  int size;

  double* samples;
  double* params;

} BandpassFilter;

void init_bpf(BandpassFilter* filter, double flow, double fhigh, int size){
  filter->lambda = M_PI * flow / (44100/2);
  filter->phi = M_PI * fhigh / (44100/2);
  filter->size = size;

  filter->samples = malloc(sizeof(double) * size);
  filter->params = malloc(sizeof(double) * size);

  for(int i = 0; i < size; i++){
    filter->samples[i] = 0;

    double dist = i - (size - 1.0) / 2.0;
    if( dist == 0.0 ){
      filter->params[i] = (filter->phi - filter->lambda) / M_PI;
    } else {
      filter->params[i] = ( sin( dist * filter->phi ) - sin( dist * filter->lambda )) / (dist * M_PI);
    }
  }
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



typedef struct {} Macro;

typedef struct {
	clock_t attack;
	clock_t release;
	double freq;
	double step;
	double phase;

	bool active;
} Note;

double bad_adsr(Note *n, clock_t now){
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
// R -- RELEASE -- cza

typedef struct {
  double a;
  double d;
  double s;
  double r;
} adsr_params;

adsr_params main_adsr_params = {0.03, 0.5, 0.2, 0.5};

double lin_adsr(Note *n, clock_t now){
  double x = (now - n->attack)/((double)CLOCKS_PER_SEC);
  double a = main_adsr_params.a, d = main_adsr_params.d,
	 s = main_adsr_params.s, r = main_adsr_params.s;
  if( n->release == -1 ){

    if( x > a + d ){
      // 1. klawisz wciśnięty długo ( więcej niż A+D ) -> etap sustain
      printf("1");
      return s;
    } else if( x < a ) {
      // 2. trzymany krócej niż A
      printf("2");
      return x / a;
    } else {
      // 3. trzymany dłużej niż A, ale krócej niż A+D
      printf("3");
      return s+(1.0-s)*(x-a)/d;
    }
    
  } else {
    if( x > a+d ){
      // klawisz puszczony po długim czasie
      // a+d-x == czas w etapie release * -1
      // /r -- do 0-1
      // *s -- do 0-s
      
      // przy okazji oznaczamy jako nieaktywny czasami
      if( (x-a-d) > r ){
	n->active = false;
	return 0.0;
      }


      return (a+d-x)*s/r;
    } else {
      // klawisz puszczony przed czasem a+d!
      // aby było gładko itd należy przez czas r wygłuszać dźwięk od
      // poprzedniej wartości do 0
      // poprzednią wartość obliczamy jakby n->release było równe -1
      // a naszym nowym x jest czas od released to tearz
      double x0 = (n->release - n->attack)/CLOCKS_PER_SEC;
      double x2 = (now - n->release)/CLOCKS_PER_SEC;
      double y;

      if( x0 < a ) {
	y = x0 / a;
      } else {
	y = s+(1.0-s)*(x0-a)/d;
      }

      return y * x2/r;

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

double (*main_adsr)(Note* note, clock_t now);
BandpassFilter main_filter;
Note all_notes[104];

static char *device = "plughw:0,0";                     /* playback device */
static snd_pcm_format_t format = SND_PCM_FORMAT_S16;    /* sample format */
static unsigned int rate = 44100;                       /* stream rate */
static unsigned int channels = 1;                       /* count of channels */
static unsigned int buffer_time = 50000;               /* ring buffer length in us */
static unsigned int period_time = 10000;               /* period time in us */
static double freq = 440;                               /* sinusoidal wave frequency in Hz */
static int verbose = 0;                                 /* verbose flag */
static int resample = 1;                                /* enable alsa-lib resampling */
static int period_event = 0;                            /* produce poll event after each period */
static snd_pcm_sframes_t buffer_size;
static snd_pcm_sframes_t period_size;
static snd_output_t *output = NULL;

static double max_phase = 2. * M_PI;

void control_loop(){

  for(int i = 0; i < 104; i++){
    all_notes[i].active = false;
    double freq = freq_calc(i);
    all_notes[i].freq = freq;
    all_notes[i].step = max_phase*(freq)/(double)rate;
  }


  struct input_event ev[1];
  int fd, rd, value, code, size = sizeof (struct input_event);
  char name[256] = "Unknown";
  char *device = NULL;

  if ((getuid ()) != 0)
    printf ("You are not root! This may not work...n/");

  device = "/dev/input/event0";
  //Open Device
  if ((fd = open (device, O_RDONLY)) == -1)
    printf ("%s is not a vaild device.n", device);
  ioctl (fd, EVIOCGNAME (sizeof (name)), name);
  printf ("IO: Reading From : %s (%s)\n", device, name);
  fflush(stdout); 


  while (1){
    if ((rd = read (fd, ev, size )) < size)
      perror("Error reading");  
    value = ev[0].value;
    code = ev[0].code;
    if (value !=2 && ev[0].type == 1){
      //got char
      printf ("Code[%d] %d \n", code, value);
      fflush(stdout);
      printf ("Code = %d \n",code);

      /*
      if( code == 12 ){
	CUTOFF -= 50;
	printf("cutoff: %lf\n", CUTOFF);
	fflush(stdout);
	continue;
      } else if (code == 13){
	CUTOFF += 50;
	printf("cutoff: %lf\n", CUTOFF);
	fflush(stdout);
	continue;
      }
      */

      int note = piano_keys[code];
      if(value == 1){
	//PRESS
	all_notes[note].phase = 0.0;
	all_notes[note].active = true;
	all_notes[note].attack = clock();
	all_notes[note].release = -1;
      }else{
	//RELEASE
	all_notes[note].release = clock();
      }
    }
  }  
  exit(EXIT_SUCCESS);
}

static void combine_sounds(const snd_pcm_channel_area_t *areas,
		snd_pcm_uframes_t offset,
		int count)
{

  unsigned char *samples[channels];
  int steps[channels];
  unsigned int chn;
  int format_bits = snd_pcm_format_width(format);
  unsigned int maxval = (1 << (format_bits - 1)) - 1;
  int bps = format_bits / 8;  /* bytes per sample */
  int phys_bps = snd_pcm_format_physical_width(format) / 8;
  int big_endian = snd_pcm_format_big_endian(format) == 1;
  int to_unsigned = snd_pcm_format_unsigned(format) == 1;
  int is_float = (format == SND_PCM_FORMAT_FLOAT_LE ||
      format == SND_PCM_FORMAT_FLOAT_BE);

  /* verify and prepare the contents of areas */
  for (chn = 0; chn < channels; chn++) {
    if ((areas[chn].first % 8) != 0) {
      printf("areas[%i].first == %i, aborting...\n", chn, areas[chn].first);
      exit(EXIT_FAILURE);
    }
    samples[chn] = /*(signed short *)*/(((unsigned char *)areas[chn].addr) + (areas[chn].first / 8));
    if ((areas[chn].step % 16) != 0) {
      printf("areas[%i].step == %i, aborting...\n", chn, areas[chn].step);
      exit(EXIT_FAILURE);
    }
    steps[chn] = areas[chn].step / 8;
    samples[chn] += offset * steps[chn];
  }


  //calc_filter_parameters();

  /* fill the channel areas */
  while (count-- > 0) {

    union {
      float f;
      int i;
    } fval;

    int res, i;
    if (is_float) {
      fval.f = sin(0);
      res = fval.i;
    } else {
      res = 0;

      clock_t clock_now = clock();
      for(int i = 28; i <= 68; i++){
	Note* n = &(all_notes[i]);

	if(! n->active )
	  continue;

	double adsr_val = main_adsr(n, clock_now);
	

	// suma z dwóch fal:
	double osc_total = osc_wave1(n->phase) * weight[0] + osc_wave2(n->phase) * weight[1];

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
    }

    if (to_unsigned)
      res ^= 1U << (format_bits - 1);

    for (chn = 0; chn < channels; chn++) {
      /* Generate data in native endian format */
      if (big_endian) {
	for (i = 0; i < bps; i++)
	  *(samples[chn] + phys_bps - 1 - i) = (res >> i * 8) & 0xff;
      } else {
	for (i = 0; i < bps; i++)
	  *(samples[chn] + i) = (res >>  i * 8) & 0xff;
      }
      samples[chn] += steps[chn];
    }
  }
}
static int set_hwparams(snd_pcm_t *handle,
		snd_pcm_hw_params_t *params,
		snd_pcm_access_t access)
{
	unsigned int rrate;
	snd_pcm_uframes_t size;
	int err, dir;
	/* choose all parameters */
	err = snd_pcm_hw_params_any(handle, params);
	if (err < 0) {
		printf("Broken configuration for playback: no configurations available: %s\n", snd_strerror(err));
		return err;
	}
	/* set hardware resampling */
	err = snd_pcm_hw_params_set_rate_resample(handle, params, resample);
	if (err < 0) {
		printf("Resampling setup failed for playback: %s\n", snd_strerror(err));
		return err;
	}
	/* set the interleaved read/write format */
	err = snd_pcm_hw_params_set_access(handle, params, access);
	if (err < 0) {
		printf("Access type not available for playback: %s\n", snd_strerror(err));
		return err;
	}
	/* set the sample format */
	err = snd_pcm_hw_params_set_format(handle, params, format);
	if (err < 0) {
		printf("Sample format not available for playback: %s\n", snd_strerror(err));
		return err;
	}
	/* set the count of channels */
	err = snd_pcm_hw_params_set_channels(handle, params, channels);
	if (err < 0) {
		printf("Channels count (%i) not available for playbacks: %s\n", channels, snd_strerror(err));
		return err;
	}
	/* set the stream rate */
	rrate = rate;
	err = snd_pcm_hw_params_set_rate_near(handle, params, &rrate, 0);
	if (err < 0) {
		printf("Rate %iHz not available for playback: %s\n", rate, snd_strerror(err));
		return err;
	}
	if (rrate != rate) {
		printf("Rate doesn't match (requested %iHz, get %iHz)\n", rate, err);
		return -EINVAL;
	}
	/* set the buffer time */
	err = snd_pcm_hw_params_set_buffer_time_near(handle, params, &buffer_time, &dir);
	if (err < 0) {
		printf("Unable to set buffer time %i for playback: %s\n", buffer_time, snd_strerror(err));
		return err;
	}
	err = snd_pcm_hw_params_get_buffer_size(params, &size);
	if (err < 0) {
		printf("Unable to get buffer size for playback: %s\n", snd_strerror(err));
		return err;
	}
	buffer_size = size;
	/* set the period time */
	err = snd_pcm_hw_params_set_period_time_near(handle, params, &period_time, &dir);
	if (err < 0) {
		printf("Unable to set period time %i for playback: %s\n", period_time, snd_strerror(err));
		return err;
	}
	err = snd_pcm_hw_params_get_period_size(params, &size, &dir);
	if (err < 0) {
		printf("Unable to get period size for playback: %s\n", snd_strerror(err));
		return err;
	}
	period_size = size;
	/* write the parameters to device */
	err = snd_pcm_hw_params(handle, params);
	if (err < 0) {
		printf("Unable to set hw params for playback: %s\n", snd_strerror(err));
		return err;
	}
	return 0;
}
static int set_swparams(snd_pcm_t *handle, snd_pcm_sw_params_t *swparams)
{
	int err;
	/* get the current swparams */
	err = snd_pcm_sw_params_current(handle, swparams);
	if (err < 0) {
		printf("Unable to determine current swparams for playback: %s\n", snd_strerror(err));
		return err;
	}
	/* start the transfer when the buffer is almost full: */
	/* (buffer_size / avail_min) * avail_min */
	err = snd_pcm_sw_params_set_start_threshold(handle, swparams, (buffer_size / period_size) * period_size);
	if (err < 0) {
		printf("Unable to set start threshold mode for playback: %s\n", snd_strerror(err));
		return err;
	}
	/* allow the transfer when at least period_size samples can be processed */
	/* or disable this mechanism when period event is enabled (aka interrupt like style processing) */
	err = snd_pcm_sw_params_set_avail_min(handle, swparams, period_event ? buffer_size : period_size);
	if (err < 0) {
		printf("Unable to set avail min for playback: %s\n", snd_strerror(err));
		return err;
	}
	/* enable period events when requested */
	if (period_event) {
		err = snd_pcm_sw_params_set_period_event(handle, swparams, 1);
		if (err < 0) {
			printf("Unable to set period event: %s\n", snd_strerror(err));
			return err;
		}
	}
	/* write the parameters to the playback device */
	err = snd_pcm_sw_params(handle, swparams);
	if (err < 0) {
		printf("Unable to set sw params for playback: %s\n", snd_strerror(err));
		return err;
	}
	return 0;
}
/*
 *   Underrun and suspend recovery
 */

static int xrun_recovery(snd_pcm_t *handle, int err)
{
	if (verbose)
		printf("stream recovery\n");
	if (err == -EPIPE) {    /* under-run */
		err = snd_pcm_prepare(handle);
		if (err < 0)
			printf("Can't recovery from underrun, prepare failed: %s\n", snd_strerror(err));
		return 0;
	} else if (err == -ESTRPIPE) {
		while ((err = snd_pcm_resume(handle)) == -EAGAIN)
			sleep(1);       /* wait until the suspend flag is released */
		if (err < 0) {
			err = snd_pcm_prepare(handle);
			if (err < 0)
				printf("Can't recovery from suspend, prepare failed: %s\n", snd_strerror(err));
		}
		return 0;
	}
	return err;
}


struct async_private_data {
	signed short *samples;
	snd_pcm_channel_area_t *areas;
	double phase;
};

clock_t last, now;

static void async_callback(snd_async_handler_t *ahandler)
{
	snd_pcm_t *handle = snd_async_handler_get_pcm(ahandler);
	struct async_private_data *data = snd_async_handler_get_callback_private(ahandler);
	signed short *samples = data->samples;
	snd_pcm_channel_area_t *areas = data->areas;
	snd_pcm_sframes_t avail;
	int err;

	last = now;
	now = clock();

	avail = snd_pcm_avail_update(handle);
	while (avail >= period_size) {
		combine_sounds(areas, 0, period_size);
		err = snd_pcm_writei(handle, samples, period_size);
		if (err < 0) {
			printf("Write error: %s\n", snd_strerror(err));
			exit(EXIT_FAILURE);
		}
		if (err != period_size) {
			printf("Write error: written %i expected %li\n", err, period_size);
			exit(EXIT_FAILURE);
		}
		avail = snd_pcm_avail_update(handle);
	}
}

static int async_loop(snd_pcm_t *handle,
		signed short *samples,
		snd_pcm_channel_area_t *areas)
{
	struct async_private_data data;
	snd_async_handler_t *ahandler;
	int err, count;
	data.samples = samples;
	data.areas = areas;
	err = snd_async_add_pcm_handler(&ahandler, handle, async_callback, &data);
	if (err < 0) {
		printf("Unable to register async handler\n");
		exit(EXIT_FAILURE);
	}
	for (count = 0; count < 2; count++) {
		combine_sounds(areas, 0, period_size);
		err = snd_pcm_writei(handle, samples, period_size);
		if (err < 0) {
			printf("Initial write error: %s\n", snd_strerror(err));
			exit(EXIT_FAILURE);
		}
		if (err != period_size) {
			printf("Initial write error: written %i expected %li\n", err, period_size);
			exit(EXIT_FAILURE);
		}
	}
	if (snd_pcm_state(handle) == SND_PCM_STATE_PREPARED) {
		err = snd_pcm_start(handle);
		if (err < 0) {
			printf("Start error: %s\n", snd_strerror(err));
			exit(EXIT_FAILURE);
		}
	}

	/* because all other work is done in the signal handler,
		 suspend the process */

	while(1) control_loop();

}
/*
 *   Transfer method - asynchronous notification + direct write
 */
static void async_direct_callback(snd_async_handler_t *ahandler)
{
	snd_pcm_t *handle = snd_async_handler_get_pcm(ahandler);
	struct async_private_data *data = snd_async_handler_get_callback_private(ahandler);
	const snd_pcm_channel_area_t *my_areas;
	snd_pcm_uframes_t offset, frames, size;
	snd_pcm_sframes_t avail, commitres;
	snd_pcm_state_t state;
	int first = 0, err;

	while (1) {
		state = snd_pcm_state(handle);
		if (state == SND_PCM_STATE_XRUN) {
			err = xrun_recovery(handle, -EPIPE);
			if (err < 0) {
				printf("XRUN recovery failed: %s\n", snd_strerror(err));
				exit(EXIT_FAILURE);
			}
			first = 1;
		} else if (state == SND_PCM_STATE_SUSPENDED) {
			err = xrun_recovery(handle, -ESTRPIPE);
			if (err < 0) {
				printf("SUSPEND recovery failed: %s\n", snd_strerror(err));
				exit(EXIT_FAILURE);
			}
		}
		avail = snd_pcm_avail_update(handle);
		if (avail < 0) {
			err = xrun_recovery(handle, avail);
			if (err < 0) {
				printf("avail update failed: %s\n", snd_strerror(err));
				exit(EXIT_FAILURE);
			}
			first = 1;
			continue;
		}
		if (avail < period_size) {
			if (first) {
				first = 0;
				err = snd_pcm_start(handle);
				if (err < 0) {
					printf("Start error: %s\n", snd_strerror(err));
					exit(EXIT_FAILURE);
				}
			} else {
				break;
			}
			continue;
		}
		size = period_size;
		while (size > 0) {
			frames = size;
			err = snd_pcm_mmap_begin(handle, &my_areas, &offset, &frames);
			if (err < 0) {
				if ((err = xrun_recovery(handle, err)) < 0) {
					printf("MMAP begin avail error: %s\n", snd_strerror(err));
					exit(EXIT_FAILURE);
				}
				first = 1;
			}
			combine_sounds(my_areas, offset, frames);
			commitres = snd_pcm_mmap_commit(handle, offset, frames);
			if (commitres < 0 || (snd_pcm_uframes_t)commitres != frames) {
				if ((err = xrun_recovery(handle, commitres >= 0 ? -EPIPE : commitres)) < 0) {
					printf("MMAP commit error: %s\n", snd_strerror(err));
					exit(EXIT_FAILURE);
				}
				first = 1;
			}
			size -= frames;
		}
	}
}
static int async_direct_loop(snd_pcm_t *handle,
		signed short *samples ATTRIBUTE_UNUSED,
		snd_pcm_channel_area_t *areas ATTRIBUTE_UNUSED)
{
	struct async_private_data data;
	snd_async_handler_t *ahandler;
	const snd_pcm_channel_area_t *my_areas;
	snd_pcm_uframes_t offset, frames, size;
	snd_pcm_sframes_t commitres;
	int err, count;
	data.samples = NULL;    /* we do not require the global sample area for direct write */
	data.areas = NULL;      /* we do not require the global areas for direct write */
	data.phase = 0;
	err = snd_async_add_pcm_handler(&ahandler, handle, async_direct_callback, &data);
	if (err < 0) {
		printf("Unable to register async handler\n");
		exit(EXIT_FAILURE);
	}
	for (count = 0; count < 2; count++) {
		size = period_size;
		while (size > 0) {
			frames = size;
			err = snd_pcm_mmap_begin(handle, &my_areas, &offset, &frames);
			if (err < 0) {
				if ((err = xrun_recovery(handle, err)) < 0) {
					printf("MMAP begin avail error: %s\n", snd_strerror(err));
					exit(EXIT_FAILURE);
				}
			}
			combine_sounds(my_areas, offset, frames);
			commitres = snd_pcm_mmap_commit(handle, offset, frames);
			if (commitres < 0 || (snd_pcm_uframes_t)commitres != frames) {
				if ((err = xrun_recovery(handle, commitres >= 0 ? -EPIPE : commitres)) < 0) {
					printf("MMAP commit error: %s\n", snd_strerror(err));
					exit(EXIT_FAILURE);
				}
			}
			size -= frames;
		}
	}
	err = snd_pcm_start(handle);
	if (err < 0) {
		printf("Start error: %s\n", snd_strerror(err));
		exit(EXIT_FAILURE);
	}
	/* because all other work is done in the signal handler,
		 suspend the process */
	while(1) control_loop();
}

/*
 *
 */
struct transfer_method {
	const char *name;
	snd_pcm_access_t access;
	int (*transfer_loop)(snd_pcm_t *handle,
			signed short *samples,
			snd_pcm_channel_area_t *areas);
};
static struct transfer_method transfer_methods[] = {
	{ "async", SND_PCM_ACCESS_RW_INTERLEAVED, async_loop },
	{ "async_direct", SND_PCM_ACCESS_MMAP_INTERLEAVED, async_direct_loop },
	{ NULL, SND_PCM_ACCESS_RW_INTERLEAVED, NULL }
};
static void help(void)
{
	int k;
	printf(
			"Usage: pcm [OPTION]... [FILE]...\n"
			"-h,--help      help\n"
			"-D,--device    playback device\n"
			"-r,--rate      stream rate in Hz\n"
			"-c,--channels  count of channels in stream\n"
			"-f,--frequency sine wave frequency in Hz\n"
			"-b,--buffer    ring buffer size in us\n"
			"-p,--period    period size in us\n"
			"-m,--method    transfer method\n"
			"-o,--format    sample format\n"
			"-v,--verbose   show the PCM setup parameters\n"
			"-n,--noresample  do not resample\n"
			"-e,--pevent    enable poll event after each period\n"
			"\n");
	printf("Recognized sample formats are:");
	for (k = 0; k < SND_PCM_FORMAT_LAST; ++k) {
		const char *s = snd_pcm_format_name(k);
		if (s)
			printf(" %s", s);
	}
	printf("\n");
	printf("Recognized transfer methods are:");
	for (k = 0; transfer_methods[k].name; k++)
		printf(" %s", transfer_methods[k].name);
	printf("\n");
}

int main(int argc, char *argv[]) {

  main_adsr = bad_adsr;
  init_bpf( &(main_filter), 3000, 4000, 17);
 

  rot_state = mmap(NULL, sizeof *rot_state, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);

  *rot_state = 100;

#ifdef MAKE_PI
  int child;

  if((child = fork()) == 0){
    prctl(PR_SET_PDEATHSIG, SIGHUP);
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

    gpioSetAlertFunc(pa, rot_callback);
    gpioSetAlertFunc(pb, rot_callback);

    while(*rot_state < 10000){sleep(10);}

    gpioSetAlertFunc(pa, 0);
    gpioSetAlertFunc(pb, 0);

    gpioTerminate();
    exit(EXIT_SUCCESS);
  }
#endif


  initscr();

  int top_row[] = { 16, 3,  17,  4, 18, 19,  6, 20,  7, 21,  8, 22 };
  int bot_row[] = { 44, 31, 45, 32, 46, 47, 34, 48, 35, 49, 36, 50 };
  init_piano_keys(52, top_row, sizeof(top_row)/sizeof(top_row[0]));
  init_piano_keys(40, bot_row, sizeof(bot_row)/sizeof(bot_row[0]));

  struct option long_option[] =
  {
    {"help", 0, NULL, 'h'},
    {"device", 1, NULL, 'D'},
    {"rate", 1, NULL, 'r'},
    {"channels", 1, NULL, 'c'},
    {"frequency", 1, NULL, 'f'},
    {"buffer", 1, NULL, 'b'},
    {"period", 1, NULL, 'p'},
    {"method", 1, NULL, 'm'},
    {"format", 1, NULL, 'o'},
    {"verbose", 1, NULL, 'v'},
    {"noresample", 1, NULL, 'n'},
    {"pevent", 1, NULL, 'e'},
    {NULL, 0, NULL, 0},
  };
  snd_pcm_t *handle;
  int err, morehelp;
  snd_pcm_hw_params_t *hwparams;
  snd_pcm_sw_params_t *swparams;
  int method = 0;
  signed short *samples;
  unsigned int chn;
  snd_pcm_channel_area_t *areas;
  snd_pcm_hw_params_alloca(&hwparams);
  snd_pcm_sw_params_alloca(&swparams);
  morehelp = 0;
  while (1) {
    int c;
    if ((c = getopt_long(argc, argv, "hD:r:c:f:b:p:m:o:vne", long_option, NULL)) < 0)
      break;
    switch (c) {
      case 'h':
        morehelp++;
        break;
      case 'D':
        device = strdup(optarg);
        break;
      case 'r':
        rate = atoi(optarg);
        rate = rate < 4000 ? 4000 : rate;
        rate = rate > 196000 ? 196000 : rate;
        break;
      case 'c':
        channels = atoi(optarg);
        channels = channels < 1 ? 1 : channels;
        channels = channels > 1024 ? 1024 : channels;
        break;
      case 'f':
        freq = atoi(optarg);
        freq = freq < 50 ? 50 : freq;
        freq = freq > 5000 ? 5000 : freq;
        break;
      case 'b':
        buffer_time = atoi(optarg);
        buffer_time = buffer_time < 1000 ? 1000 : buffer_time;
        buffer_time = buffer_time > 1000000 ? 1000000 : buffer_time;
        break;
      case 'p':
        period_time = atoi(optarg);
        period_time = period_time < 1000 ? 1000 : period_time;
        period_time = period_time > 1000000 ? 1000000 : period_time;
        break;
      case 'm':
        for (method = 0; transfer_methods[method].name; method++)
          if (!strcasecmp(transfer_methods[method].name, optarg))
            break;
        if (transfer_methods[method].name == NULL)
          method = 0;
        break;
      case 'o':
        for (format = 0; format < SND_PCM_FORMAT_LAST; format++) {
          const char *format_name = snd_pcm_format_name(format);
          if (format_name)
            if (!strcasecmp(format_name, optarg))
              break;
        }
        if (format == SND_PCM_FORMAT_LAST)
          format = SND_PCM_FORMAT_S16;
        if (!snd_pcm_format_linear(format) &&
            !(format == SND_PCM_FORMAT_FLOAT_LE ||
              format == SND_PCM_FORMAT_FLOAT_BE)) {
          printf("Invalid (non-linear/float) format %s\n",
              optarg);
          return 1;
        }
        break;
      case 'v':
        verbose = 1;
        break;
      case 'n':
        resample = 0;
        break;
      case 'e':
        period_event = 1;
        break;
    }
  }
  if (morehelp) {
    help();
    return 0;
  }
  err = snd_output_stdio_attach(&output, stdout, 0);
  if (err < 0) {
    printf("Output failed: %s\n", snd_strerror(err));
    return 0;
  }
  printf("Playback device is %s\n", device);
  printf("Stream parameters are %iHz, %s, %i channels\n", rate, snd_pcm_format_name(format), channels);
  printf("Sine wave rate is %.4fHz\n", freq);
  printf("Using transfer method: %s\n", transfer_methods[method].name);
  if ((err = snd_pcm_open(&handle, device, SND_PCM_STREAM_PLAYBACK, 0)) < 0) {
    printf("Playback open error: %s\n", snd_strerror(err));
    return 0;
  }

  if ((err = set_hwparams(handle, hwparams, transfer_methods[method].access)) < 0) {
    printf("Setting of hwparams failed: %s\n", snd_strerror(err));
    exit(EXIT_FAILURE);
  }
  if ((err = set_swparams(handle, swparams)) < 0) {
    printf("Setting of swparams failed: %s\n", snd_strerror(err));
    exit(EXIT_FAILURE);
  }
  if (verbose > 0)
    snd_pcm_dump(handle, output);
  samples = malloc((period_size * channels * snd_pcm_format_physical_width(format)) / 8);
  if (samples == NULL) {
    printf("No enough memory\n");
    exit(EXIT_FAILURE);
  }

  areas = calloc(channels, sizeof(snd_pcm_channel_area_t));
  if (areas == NULL) {
    printf("No enough memory\n");
    exit(EXIT_FAILURE);
  }
  for (chn = 0; chn < channels; chn++) {
    areas[chn].addr = samples;
    areas[chn].first = chn * snd_pcm_format_physical_width(format);
    areas[chn].step = channels * snd_pcm_format_physical_width(format);
  }

  err = transfer_methods[method].transfer_loop(handle, samples, areas);

  if (err < 0)
    printf("Transfer failed: %s\n", snd_strerror(err));
  free(areas);
  free(samples);
  snd_pcm_close(handle);

#ifdef MAKE_PI
  kill(child, 9);
#endif
  wait(NULL);
  munmap(rot_state, sizeof *rot_state);

  endwin();
  return 0;
}
