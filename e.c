/*
 *  This small demo sends a simple sinusoidal wave to your speakers.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sched.h>
#include <errno.h>
#include <getopt.h>
#include <alsa/asoundlib.h>
#include <sys/time.h>
#include <math.h>

#include <ncurses.h>



typedef struct {} Macro;

typedef struct {
  float (*adsr)(Note);
} Sound;

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
  double seconds = 0.33;
  double r = fmin(1.0, fmax(0.0, 1.0 - ((double)(now - n->release))/CLOCKS_PER_SEC/(seconds)));
  if ( r == 0.0 )
    n->active = false;
  return r;
}

double freq_calc(int n){ // A4 = 49 -> 440Hz
  return 440.0 * pow(2.0, (n-49)/12.0);
}

Note all_notes[104];

Sound main_sound = { bad_adsr };

static char *device = "plughw:00,00";                     /* playback device */
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

  while (1) {
    char c;
    scanf("%c\n",&c);
    printf("bing %c\n", c);
    while( c != 'x' ){
      printf("%i\n", c);
      scanf("%c\n", &c);

      if( c >= 'a' && c <= 'g' ){
        int note = 49 + c - 'a';
				int used = 0;
				for(int i = 40; i < 49+13; i++){
								if(all_notes[i].active)
												used += 1;
				}
				printf("used: %d\n", used);
        printf("Note: %i\n", note);

				all_notes[note].phase = 0.0;
				all_notes[note].active = true;
				all_notes[note].attack = clock();
				all_notes[note].release = clock()+0.08*CLOCKS_PER_SEC;

				all_notes[note+5].phase = 0.0;
				all_notes[note+5].active = true;
				all_notes[note+5].attack = clock();
				all_notes[note+5].release = clock()+0.08*CLOCKS_PER_SEC;

				all_notes[note+7].phase = 0.0;
				all_notes[note+7].active = true;
				all_notes[note+7].attack = clock();
				all_notes[note+7].release = clock()+0.08*CLOCKS_PER_SEC;


      }
    }
    exit(EXIT_SUCCESS);
  }

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

  double whatwave(double x){
    if(x > M_PI)
						return 1.0;
		else
						return -1.0;
  }

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
      printf("IS_FLOAT!!\n");
    } else {
      res = 0;

			clock_t clock_now = clock();
			for(int i = 49-13; i < 49+26; i++){
				Note* n = &(all_notes[i]);
				if(! n->active )
					continue;
				double adsr_val = bad_adsr(n, clock_now);
				res += whatwave(n->phase) * maxval * adsr_val;

        n->phase += n->step;
        if (n->phase >= max_phase){
          n->phase -= max_phase;
				}
      }

			res /= 20.0;
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

int main(int argc, char *argv[])
{

  //initscr();

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

  //endwin();
  return 0;
}

