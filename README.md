# vinst2sfz
Record MIDI instrument and produce SFZ mapping and WAV samples
-----------------------------------
# Description
Record audio samples from virtual instrument (software or hardware) and generate SFZ mapping for them  
		-> PCM audio data are taken from specified JACK output  
		-> the instrument is triggered through ALSA MIDI sequencer  
		-> it is intended primary for drum samples    
                -> Linux only  

# Usage
vinst2sfz -m file-name-base -n notenumber -p jack-port-string [options]  

   -a target ALSA MIDI port (Midi-Through 14:0)  
   -c target ALSA MIDI channel (10=Drums)  
   -d maximal duration in seconds (10 s)  
   -h help  
   -i input threshold instant value (-50 dB)  
   -l list ALSA MIDI ports  
   -o output threshold RMS value(-70 dB)  
   -r number of round robins (1)  
   -t target dynamic range (20 dB)  
   -v number of velocity layers (5)  
   -x dynamic range and velocity layer boundary analysis only -> no file output  
   -y velocity curve correction in SFZ  
      ----> 0=no correction (default)  
      ----> 1=under specified target range normalized to target range and above without correction  
      ----> 2=original dynamic range divided to 'v' equal steps  
      ----> 3=specified target dynamic range divided to 'v' equal steps  
  
# Build
gcc -lm -lasound -lsndfile -ljack -pthread -o vinst2sfz vinst2sfz.c

# Author
Miroslav Kovac (mixxoo@gmail.com)
