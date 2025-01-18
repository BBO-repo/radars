# radars
Illustration of some radar concepts with octave/matlab:
1. radar_1_pulse_radar.m simulate pulse radar with matched filter processing and CFAR detection to retrieve target position.
2. radar_2_pulse_integration.m similar to 1. with additional non coherent pulse integration.
3. radar_3_pulse_compression.m illustrates pulse compression with linear frequency modulation (LFM).
4. radar_3b_pulse_compression.m illustrates range resolution gain of rectangular pulse versus LFM pulse compression - two close targets are detected with LFM while not with rectangular pulse for the same pulse duration.
5. radar_4_continous_wave_LFM illustrates a range doppler map of a moving target with 2D FFT. 

# requirements
If run from matlab, make sure you have the signal processing toolbox.<br>
If run from octave, make sure signal package is installer.<br>
If using the joined Dockerfile, it will work out of the box.<br>

# how to use
From an matlab/octave command prompt just run the desired script, for example:
```
>> radar_1_pulse_radar
```
If using docker, making X11 apps work inside docker by running this command in host
```
xhost +local:docker
```
Then inside the docker container running the following command allow to open the octave gui from which you can run any script:
```
octave --gui
```
Inside the matlab/octave gui command prompt, you can run:
```
>> radar_1_pulse_radar
```