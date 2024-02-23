/********************************************************************************************************
	Record audio samples from virtual instrument (software or hardware) and generate SFZ mapping for them
		-> PCM audio data are taken from specified JACK output
		-> the instrument is triggered through ALSA MIDI sequencer
		-> it is intended primary for drum samples

		Author : Miroslav Kovac (mixxoo@gmail.com)
*********************************************************************************************************/
#include <signal.h>
#include <pthread.h>
#include <math.h>
#include <alsa/asoundlib.h>
#include <jack/jack.h>
#include <sndfile.h>

static snd_seq_t *seq;

static snd_seq_event_t ev;
static char *dname;			/* file name base for naming WAVs and SFZ file*/
static char *jports;			/* JACK output port string to connect to */
static char *addr_s;			/* ALSA SEQ address string to connect to */
static uint8_t note = 36;	/* MIDI note number */
static uint8_t vel = 5;		/* number of velocity layers */
static uint8_t rrob = 1;	/* number of round-robins*/
static uint8_t midich = 9; //MIDI channel zero based -> ch 10=Drums
static uint8_t dur = 10;	//seconds -> sets how much memory we need to allocate
static float target_d_range = 20.0f; /* target dynamic range in dB */

/* Amplitude correction type :
   ----> 0=no correction : SFZ should reproduce the original velocity->amplitude mapping
   ----> 1=amplitudes under target range are normalized to target range and
           amplitudes above target range are left without correction 
   ----> 2=measured dynamic range is divided into equal steps according to number of velocity layers
   ----> 3=target dynamic range is divided into equal steps according to number of velocity layers*/
static char ampl_correction = 0;

static const size_t sample_size = sizeof(jack_default_audio_sample_t);
static jack_client_t *client;
static jack_port_t *inp1,*inp2;
static jack_default_audio_sample_t *ins1,*ins2;		/* JACK input buffers */
static jack_default_audio_sample_t *smpl1,*smpl2;	/* sample memory pointers */
static jack_default_audio_sample_t *samples;			/* current sample memory pointer */

static double inthr = 0.003162277660168;				/* -50 dB */
static double outthr = 0.000316227766016;				/* -70 dB */
static jack_nframes_t fsize, wrt;						/* number of allocated and written samples */
static jack_nframes_t rsmpln;								/* number of read samples */
static jack_nframes_t srate;								/* JACK server sample rate */
static jack_nframes_t bufsize;							/* JACK buffer size */

static float *s_max_vals;									/* max values of recorded samples */

static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t data_ready = PTHREAD_COND_INITIALIZER;

static char wavname[256];				/* output WAV file name */
static char sfzname[256];				/* output SFZ file name */

/**
 * The process callback for this JACK application is called in a
 * special realtime thread once for each audio cycle.
 *
 */
int process (jack_nframes_t nframes, void *arg)
{
	static char ready = 1;				/* indicates whether the process thread is ready */
	static double rms_sum1 = 0.0;		/* moving RMS value channel 1*/
	static double rms_sum2 = 0.0;		/* moving RMS value channel 2*/
	static jack_nframes_t rms_size;	/* window size for moving RMS */
	static jack_nframes_t rms_start;	/* start of the window for moving RMS */
	static jack_nframes_t rstart;		/* whether recording has started */

	if (ready) {
		if (pthread_mutex_trylock(&lock) == 0) {
			wrt = 0;
			rsmpln = 0;
			ready = 0;
			rstart = 0;
			rms_sum1 = 0.0;
			rms_sum2 = 0.0;
			rms_size = srate / 20; /* 50 ms */
			rms_start = 0;
		} else return 0;
	}

	if (wrt < fsize) {
		/* record until there is memory available */
		ins1 = jack_port_get_buffer(inp1, nframes);
		ins2 = jack_port_get_buffer(inp2, nframes);

		jack_nframes_t i;

		for (i=0; i < nframes;++i) {

			/* read sample counter */
			rsmpln++;
			rsmpln++;

			if (wrt == 0) {
				/* Recording starts with the first sample above input threshold */
				if (fabsf(ins1[i]) > inthr) ++rstart;
				if (fabsf(ins2[i]) > inthr) ++rstart;
			}
			if (rstart > 0) {
				/* store samples into memory -> 2 channels interleaved */
				samples[wrt++] = ins1[i];
				samples[wrt++] = ins2[i];
				/* add to moving RMS sums */
				rms_sum1 += ins1[i] * ins1[i];
				rms_sum2 += ins2[i] * ins2[i];
				/* subtract from moving RMS sums */
				if (wrt > rms_size * 2) {
					rms_sum1 -= samples[rms_start] * samples[rms_start++];
					rms_sum2 -= samples[rms_start] * samples[rms_start++];
				}
			}
		}

		if (rsmpln > rms_size) {

			/* Calculate RMS values */

			double rms_val1 = sqrt(rms_sum1 / rms_size);
			double rms_val2 = sqrt(rms_sum2 / rms_size);

			/* Until both channels RMS values are above output threshold we do not finish recording */
			if (rms_val1 > outthr || rms_val2 > outthr) return 0;

		} else return 0;
	}
	/* End of recording a single strike */
	ready = 1;
	pthread_cond_signal(&data_ready);
	pthread_mutex_unlock(&lock);
	return 0;
}

/**
 * JACK calls this shutdown_callback if the server ever shuts down or
 * decides to disconnect the client.
 */
void jack_shutdown (void *arg)
{
	exit (EXIT_FAILURE);
}

static void list_ports(void)
{

	int err;

	/* open sequencer */
	err = snd_seq_open(&seq, "default", SND_SEQ_OPEN_DUPLEX, 0);

	snd_seq_client_info_t *cinfo;
	snd_seq_port_info_t *pinfo;

	snd_seq_client_info_alloca(&cinfo);
	snd_seq_port_info_alloca(&pinfo);

	puts(" Port    Client name                      Port name");

	snd_seq_client_info_set_client(cinfo, -1);
	while (snd_seq_query_next_client(seq, cinfo) >= 0) {
		int client = snd_seq_client_info_get_client(cinfo);

		snd_seq_port_info_set_client(pinfo, client);
		snd_seq_port_info_set_port(pinfo, -1);
		while (snd_seq_query_next_port(seq, pinfo) >= 0) {
			/* port must understand MIDI messages */
		/*	if (!(snd_seq_port_info_get_type(pinfo)
		     & SND_SEQ_PORT_TYPE_MIDI_GENERIC))
		continue;
			/* we need both WRITE and SUBS_WRITE */
			/*if ((snd_seq_port_info_get_capability(pinfo)
		    & (SND_SEQ_PORT_CAP_WRITE | SND_SEQ_PORT_CAP_SUBS_WRITE))
		   != (SND_SEQ_PORT_CAP_WRITE | SND_SEQ_PORT_CAP_SUBS_WRITE))
		continue;*/
			printf("%3d:%-3d  %-32.32s %s\n",
		      snd_seq_port_info_get_client(pinfo),
		      snd_seq_port_info_get_port(pinfo),
		      snd_seq_client_info_get_name(cinfo),
		      snd_seq_port_info_get_name(pinfo));
		}
	}

	snd_seq_close(seq);
}

void signal_handler(int sig)
{
	fprintf(stderr, "\nsignal received, exiting ...\n");
	exit(EXIT_SUCCESS);
}

void init_alsa()
{
	/*************************************/
	/* Create an connect ALSA SEQ client */
	/*************************************/
	int err = snd_seq_open(&seq, "default", SND_SEQ_OPEN_OUTPUT, 0);
	if (err < 0) {
		printf("Error %d\n",err);
		exit (EXIT_FAILURE);
	}

	snd_seq_set_client_name(seq, "vinst2sfz");

	err = snd_seq_create_simple_port(seq, 
                                    "Port 0",
                                    SND_SEQ_PORT_CAP_READ|SND_SEQ_PORT_CAP_SUBS_READ,
                                    SND_SEQ_PORT_TYPE_MIDI_GENERIC);
	if (err < 0)
	{
		printf("Error %d\n",err);
		exit (EXIT_FAILURE);
	}
	else
	{
		printf("Created ALSA SEQ port \'vinst2sfz\' : %d\n",err);
	}

	if (addr_s == NULL) {
		//Connect to MIDI Through
		snd_seq_connect_to(seq, 0, 14, 0);
	} else {
		/* Subscribe to selected destination ALSA MIDI port */
		snd_seq_addr_t addr;
		snd_seq_parse_address(seq,&addr,addr_s);
		snd_seq_connect_to(seq, 0, addr.client, addr.port);
	}

	/* Set sequencer event defaults */

	snd_seq_ev_clear(&ev);
	snd_seq_ev_set_source(&ev, 0);
	snd_seq_ev_set_subs(&ev);
	snd_seq_ev_set_direct(&ev);
}

void init_jack()
{
	/**************************/
	/* JACK interface setting */
	/**************************/

	const char *client_name = "vinst2sfz";
	const char *server_name = NULL;
	jack_options_t options = JackNoStartServer;
	jack_status_t status;

	/* open a client connection to the JACK server */

	client = jack_client_open (client_name, options, &status, server_name);
	if (client == NULL) {
		fprintf (stderr, "jack_client_open() failed, "
			 "status = 0x%2.0x\n", status);
		if (status & JackServerFailed) {
			fprintf (stderr, "Unable to connect to JACK server\n");
		}
		exit (EXIT_FAILURE);
	}

	if (status & JackNameNotUnique) {
		client_name = jack_get_client_name(client);
		fprintf (stderr, "unique name `%s' assigned\n", client_name);
	}

	/* tell the JACK server to call `process()' whenever
	   there is work to be done.
	*/

	jack_set_process_callback (client, process, 0);

	/* tell the JACK server to call `jack_shutdown()' if
	   it ever shuts down, either entirely, or if it
	   just decides to stop calling us.
	*/

	jack_on_shutdown (client, jack_shutdown, 0);

	/* display the current sample rate and buffer size. */

	srate = jack_get_sample_rate (client);
	bufsize = jack_get_buffer_size (client);
	printf ("JACK engine sample rate: %u\n", srate);
	printf ("JACK buffer size : %u\n", bufsize);

	/* create input ports */

	inp1 = jack_port_register (client, "1",
					 JACK_DEFAULT_AUDIO_TYPE,
					 JackPortIsInput, 0);

	if (inp1 == NULL) {
		fprintf(stderr, "no more JACK ports available\n");
		exit (EXIT_FAILURE);
	}

	inp2 = jack_port_register (client, "2",
					 JACK_DEFAULT_AUDIO_TYPE,
					 JackPortIsInput, 0);

	if (inp2 == NULL) {
		fprintf(stderr, "no more JACK ports available\n");
		exit (EXIT_FAILURE);
	}

	/* lock the mutex to finish initialization */
	pthread_mutex_lock (&lock);

	/* Tell the JACK server that we are ready to roll.  Our
	 * process() callback will start running now. */

	if (jack_activate (client)) {
		fprintf (stderr, "cannot activate client");
		exit (EXIT_FAILURE);
	}

	/* Connect the ports.  You can't do this before the client is
	 * activated, because we can't make connections to clients
	 * that aren't running.  Note the confusing (but necessary)
	 * orientation of the driver backend ports: playback ports are
	 * "input" to the backend, and capture ports are "output" from
	 * it.
	 */

	const char **ports;
	ports = jack_get_ports (client, jports, NULL, JackPortIsOutput);
	if (ports == NULL) {
		fprintf(stderr, "no required JACK ports found\n");
		exit (EXIT_FAILURE);
	}

	if (jack_connect (client, ports[0], jack_port_name (inp1))) {
		fprintf (stderr, "cannot connect input port 1\n");
	}

	if (jack_connect (client, ports[1], jack_port_name (inp2))) {
		fprintf (stderr, "cannot connect input port 2\n");
	}

	free (ports);

}

void analyze()
{
	/* Samples are stored only in memory.
		Cross-correlation coefficients are computed between neighboring velocity levels
		to estimate where is the velocity level boundary in source device and so to suggest
		an appropriate number of velocity layers needed to represent original dynamic range
		as close as possible.
	*/

	/************************************************************************/
	/* Prepare required memory for recorded samples                         */
	/************************************************************************/

	/* number of frames for 2 channels and required length and sample rate */
	fsize = dur * srate * 2;

	/* number of frames must be a whole number multiplier of (channel number) x (bufsize) */
	fsize /= bufsize * 2;
	fsize *= bufsize * 2;

	size_t fsize_b = fsize * sample_size;

	smpl1 = malloc(fsize_b);
	smpl2 = malloc(fsize_b);

	samples = smpl1;

	printf("Allocated memory size : %u bytes\n",fsize_b*2);

	uint8_t vels = 127/vel;		/* velocity step */
	uint8_t avel = 0;				/* actual MIDI velocity */
	uint8_t pvel = 0;				/* previous MIDI velocity */
	jack_nframes_t pwrt = 0;	/* previously written number of samples */
	jack_nframes_t mwrt = 0;	/* minimal of written and previously written samples */

	while (avel < 127) {

		pvel = avel;
		avel += vels;
		if (avel+vels > 127) avel = 127;

		pwrt = wrt;

		/* Send NOTE ON */
		snd_seq_ev_set_noteon(&ev,midich,note,avel);
		snd_seq_event_output(seq, &ev);
		snd_seq_drain_output(seq);

		/* Wait for JACK process to record required audio data */
		pthread_cond_wait(&data_ready, &lock);

		/* Send NOTE OFF */
		snd_seq_ev_set_noteoff(&ev,midich,note,avel);
		snd_seq_event_output(seq, &ev);
		snd_seq_drain_output(seq);

		if (pvel > 0) {
			/* we have at least 2 vectors recorded so
				we can calculate correlation coefficient */

			/* get number of used samples 
				both vectors must have equal length */
			if (wrt < pwrt)
				mwrt = wrt;
			else
				mwrt = pwrt;

			/* maximal amplitudes */
			float mx1 = 0.0f;
			float mx2 = 0.0f;

			/* vector sums for computation of correlation coefficient */
			double s12 = 0.0f;
			double s11 = 0.0f;
			double s22 = 0.0f;

			for (jack_nframes_t k = 0;k<=mwrt;k++) {
				s12 += smpl1[k] * smpl2[k];
				s11 += smpl1[k] * smpl1[k];
				s22 += smpl2[k] * smpl2[k];

				jack_default_audio_sample_t act_val;

				act_val = fabsf(smpl1[k]);
				if (act_val > mx1) mx1 = act_val;

				act_val = fabsf(smpl2[k]);
				if (act_val > mx2) mx2 = act_val;
			}

			/* get maximals in dB*/
			float mx1db;
			float mx2db;

			if (samples == smpl2) {
				mx1db = 20.0*log10f(mx1);
				mx2db = 20.0*log10f(mx2);
			} else {
				mx1db = 20.0*log10f(mx2);
				mx2db = 20.0*log10f(mx1);
			}

			/* correlation coefficient */
			double xcc = s12 / sqrt(s11 * s22);

			printf("Vel %.3u <--> Vel %.3u : XCC = %.5f : Peak value = %.1f dB <--> %.1f dB\n",pvel,avel,xcc,mx1db,mx2db);

		}

		/* Exchange buffers */
		if (samples == smpl1)
			samples = smpl2;
		else
			samples = smpl1;

	}
}

void record_and_gensfz()
{
	/************************************************************************/
	/* Prepare required memory for recorded samples                         */
	/************************************************************************/

	/* number of frames for 2 channels and required length and sample rate */
	fsize = dur * srate * 2;

	/* number of frames must be a whole number multiplier of (channel number) x (bufsize) */
	fsize /= bufsize * 2;
	fsize *= bufsize * 2;

	size_t fsize_b = fsize * sample_size;

	samples = malloc(fsize_b);

	printf("Allocated file size : %u bytes\n",fsize_b);

	/************************************************************************/
	/* Set parameters of output WAV file                                    */
	/************************************************************************/
	SNDFILE *sf;
	SF_INFO sf_info;
	sf_count_t sf_wrt_b;

	sf_info.samplerate = srate;
	sf_info.channels = 2;
	sf_info.format = SF_FORMAT_WAV|SF_FORMAT_PCM_24;
	/************************************************************************/

	FILE *sfzfile;	/* Pointer to SFZ text file */

	uint8_t vels = 127/vel;		/* velocity step */
	uint8_t hivel = 0;			/* actual MIDI velocity */
	uint8_t r;						/* Round robin counter */
	float sv_min = 1.0f;			/* minimal maximum of all recording */
	float sv_max = 0.0f;			/* maximum amplitude of all recording */
	float sv_d_range;				/* dynamic range of all recording in dB */
	uint16_t scounter = vel*rrob;	/* sample counter */

	/* Allocate memory for sample maximal values */
	s_max_vals = calloc(scounter, sizeof(float));
	printf("Awaiting %u sample files\n",scounter);
	scounter = 0;

	sprintf(sfzname,"%s.sfz",dname);

	/* open SFZ and write the header */
	if ((sfzfile = fopen(sfzname,"w")) == NULL) {
		fprintf (stderr, "cannot open sfzfile for write (%s)\n", sfzname);
		exit (EXIT_FAILURE);
	}

	fputs("//***************************************************\n",sfzfile);
	fprintf(sfzfile,"//   %s\n",dname);
	fputs("//***************************************************\n",sfzfile);
	fputs("//           generated by vinst2sfz\n",sfzfile);
	fputs("//***************************************************\n",sfzfile);

	/*********************************************************************/
	/* Loop through all velocity layers and round robins and save WAVs   */
	/*********************************************************************/
	while (hivel < 127) {
		hivel += vels;
		if (hivel+vels > 127) hivel = 127;

		for (r=1;r<=rrob;r++) {

			if (rrob > 1)
				sprintf(wavname,"%s-%.3u-%.2u.wav",dname,hivel,r);
			else
				sprintf(wavname,"%s-%.3u.wav",dname,hivel);

			printf("%s",wavname);
			fflush(stdout);
			if ((sf = sf_open(wavname,SFM_WRITE,&sf_info)) == NULL) {
				char errstr[256];
				sf_error_str (0, errstr, sizeof (errstr) - 1);
				fprintf (stderr, "cannot open sndfile for output (%s)\n", errstr);
				exit (EXIT_FAILURE);
			}

			/* Send NOTE ON */
			snd_seq_ev_set_noteon(&ev,midich,note,hivel);
			snd_seq_event_output(seq, &ev);
			snd_seq_drain_output(seq);

			/* Wait for JACK process to record required audio data */
			pthread_cond_wait(&data_ready, &lock);

			/* Send NOTE OFF */
			snd_seq_ev_set_noteoff(&ev,midich,note,hivel);
			snd_seq_event_output(seq, &ev);
			snd_seq_drain_output(seq);

			/* Get dynamic range of recorded sample to be used in SFZ mapping */
			float max_sv = 0.0f;

			for (jack_nframes_t i = 0; i < wrt;++i) {
				jack_default_audio_sample_t act_val = fabsf(samples[i]);
				if (act_val > max_sv) max_sv = act_val;
			}
			/* Store for overall minimal maximum and maximal maximum */
			if (max_sv < sv_min) sv_min = max_sv;
			if (max_sv > sv_max) sv_max = max_sv;

			/* Write output WAV file */
			sf_wrt_b = sf_write_float(sf, samples, wrt);
			sf_close(sf);

			/* Convert to dB and write to the array, output and SFZ file */
			max_sv = 20.0*log10f(max_sv);
			s_max_vals[scounter++] = max_sv;
			fprintf(sfzfile,"//   %s = %.2f dB\n",wavname,max_sv);
			printf(" DONE   Peak value : %.2f dB\n",max_sv);
		}
	}
	sv_d_range = 20.0*log10f(sv_max/sv_min);	/* overall dynamic range in dB */
	printf("Overall dynamic range : %.1f dB\n",sv_d_range);

	/*******************************************************/
	/*               Generate SFZ file                     */
	/*******************************************************/
	fputs("//***************************************************",sfzfile);
	if (rrob > 1)
		fprintf(sfzfile,"\n\n<global> loop_mode=one_shot key=%u",note);
	else
		fprintf(sfzfile,"\n\n<group> loop_mode=one_shot key=%u",note);

	/* Reloop the cycle one more and fill the content of SFZ */
	uint8_t lovel = 1;
	float d_range_step;		//dynamic range step
	float d_range_target;	//dynamic range target
	float pd_range_target;	//previous dynamic range target

	if (ampl_correction == 2) {

		d_range_step = sv_d_range / vel;
		d_range_target = -sv_d_range;
		pd_range_target = d_range_target - 3.0f;

	} else if (ampl_correction == 3) {

		d_range_step = target_d_range / vel;
		d_range_target = -target_d_range;
		pd_range_target = d_range_target - 3.0f;

	}

	float velcurve = 0.707f;

	hivel = 0;
	scounter = 0;

	while (hivel < 127) {
		hivel += vels;
		if (hivel+vels > 127) hivel = 127;

		if (ampl_correction > 1) {

			if (hivel == 127) d_range_target = 0.0f;
			velcurve = powf(10,(pd_range_target-d_range_target)/20.0f);

		} else {

			d_range_target = s_max_vals[scounter];

			if (scounter > 0)
				velcurve = powf(10,(s_max_vals[scounter-1]-s_max_vals[scounter])/20.0f);

			if (ampl_correction == 1 && d_range_target < -target_d_range) {
				d_range_target = -target_d_range;
				velcurve = 1.0f;
			}

		}

		float volume;

		if (rrob > 1) {

			fprintf(sfzfile,"\n\n<group> lovel=%u hivel=%u amp_velcurve_%u=%.3f amp_velcurve_%u=1",
							lovel,hivel,lovel,velcurve,hivel);
			fprintf(sfzfile,"\n//   Target volume %.1f dB",d_range_target);

			float rand_step = 1.0f/rrob;
			float lorand = 0.0f;
			float hirand = rand_step;

			for (r=1;r<=rrob;r++) {
				volume = d_range_target - s_max_vals[scounter++];
				if (r == rrob) hirand = 1.0f;

				fprintf(sfzfile,"\n<region> sample=%s-%.3u-%.2u.wav volume=%.2f lorand=%.2f hirand=%.2f",
							dname,hivel,r,volume,lorand,hirand);

				lorand = hirand;
				hirand += rand_step;
			}

		} else {
			//single sample file per velocity layer -> group not needed
			volume = d_range_target - s_max_vals[scounter++];

			fprintf(sfzfile,"\n<region> sample=%s-%.3u.wav lovel=%u hivel=%u volume=%.2f amp_velcurve_%u=%.3f amp_velcurve_%u=1",
						dname,hivel,lovel,hivel,volume,lovel,velcurve,hivel);
			fprintf(sfzfile,"  //Target volume %.1f dB",d_range_target);
		}

		lovel = hivel+1;

		if (ampl_correction > 1) {
			pd_range_target = d_range_target;
			d_range_target += d_range_step;
		}
		
	}

	fclose(sfzfile);

}

void show_usage()
{
	puts("\nvinst2sfz -m file-name-base -n notenumber -p jack-port-string [options]\n");
	puts("   -a target ALSA MIDI port (Midi-Through 14:0)");
	puts("   -c target ALSA MIDI channel (10=Drums)");
	puts("   -d maximal duration in seconds (10 s)");
	puts("   -h help");
	puts("   -i input threshold instant value (-50 dB)");
	puts("   -l list ALSA MIDI ports");
	puts("   -o output threshold RMS value(-70 dB)");
	puts("   -r number of round robins (1)");
	puts("   -t target dynamic range (20 dB)");
	puts("   -v number of velocity layers (5)");
	puts("   -x dynamic range and velocity layer boundary analysis only -> no file output");
	puts("   -y velocity curve correction in SFZ");
	puts("      ----> 0=no correction (default)");
	puts("      ----> 1=under specified target range normalized to target range and above without correction");
	puts("      ----> 2=original dynamic range divided to \'v\' equal steps");
	puts("      ----> 3=specified target dynamic range divided to \'v\' equal steps\n");
}

void cleanup()
{
	snd_seq_close(seq);
	jack_client_close(client);
	free(samples);
	free(s_max_vals);
	puts("All done. Bye!\n");
}

int main (int argc, char *argv[])
{
	/*******************/
	/* Initialization  */
	/*******************/
	char do_list = 0;
	char do_analysis = 0;

	opterr = 0;
	int opt,err;
	char optok = 0;
	while ((opt = getopt(argc, argv, "a:d:hi:lm:n:o:p:r:t:v:xy:")) != -1) {
		switch (opt) {
			case 'a':
				addr_s = optarg;
				break;
			case 'c':
				midich = atoi(optarg)-1;
				break;
			case 'd':
			   dur = atoi(optarg);
			   break;
			case 'h':
			   show_usage();
			   exit(EXIT_SUCCESS);
			case 'i':
			   inthr = pow(10.0,-atof(optarg)/20.0);
			   break;
			case 'l':
				do_list = 1;
				break;
			case 'm':
			   dname = optarg;
				++optok;
			   break;
			case 'n':
			   note = atoi(optarg);
				++optok;
			   break;
			case 'o':
			   outthr = pow(10.0,-atof(optarg)/20.0);
			   break;
			case 'p':
			   jports = optarg;
				++optok;
			   break;
			case 'r':
			   rrob = atoi(optarg);
			   break;
			case 't':
			   target_d_range = strtof(optarg, NULL);
			   break;
			case 'v':
			   vel = atoi(optarg);
			   break;
			case 'x':
			   do_analysis = 1;
				++optok;
			   break;
			case 'y':
			   ampl_correction = atoi(optarg);
			   break;
			default: /* '?' */
			   show_usage();
			   exit(EXIT_FAILURE);
		}
	}

	if (do_list) {
		list_ports();
		exit(EXIT_SUCCESS);
	}

	if (argc < 7 || optok < 3) {
	   show_usage();
	   exit(EXIT_FAILURE);
	}

	err = atexit(cleanup);
	if (err != 0)
	{
		fprintf(stderr, "cannot set exit function\n");
        exit(EXIT_FAILURE);
	}

	/************************************************************************/
	signal(SIGTERM, signal_handler);
	signal(SIGINT, signal_handler);
	/************************************************************************/

	init_alsa();
	init_jack();

	if (do_analysis) 
		analyze();
	else
		record_and_gensfz();

	exit(EXIT_SUCCESS);
}