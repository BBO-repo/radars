clc;
clear all;
close all;

if ( exist('OCTAVE_VERSION', 'builtin') ~= 0 )
  % load signal package
  run .octaverc
endif

%% radar design parameters
max_range = 200;        % [m] maximum detection range
range_resolution = 1;   % [m] range resolution

%% generate FMCW Waveform signal
c = 3e8;                            % speed of light
chirp_duration = 6*(2*max_range/c); % [s] single chirp time duration
bandwidth = c/(2*range_resolution); % [Hz] FMCW sweep bandwidth
slope = bandwidth/chirp_duration;   % [Hz/s] slope of the linear FMCW

fc = 10e9;     % [Hz] radar carrier frequency
lambda = c/fc; % [m] radar wavelength

n_pulses = 128;    % number of sent chirps
n_samples = 1024;  % number of samples in one chirp (time vector or range cells)
fast_time = chirp_duration*((0:n_samples-1)/n_samples).';   % time vector for a single chirp
fs = chirp_duration/n_samples;                              % [Hz] sampling frequency

% build time vector [fast_time fast_time+chirp_duration fast_time+ 2*chirp_duration ..]
sim_time = ones(n_samples,1)*(0:n_pulses-1)*chirp_duration; % time array n_samples x n_pulses
sim_time = sim_time + fast_time*ones(1,n_pulses);           % a

% generate transmitted linear chirp
tx_signal = exp(1i*2*pi*(fc*fast_time + slope * (fast_time.^2)/2)); % transmitted signal

% generate received signal
init_range   = 170;   % Initial distance of the target
target_speed = 100;    % Velocity of the target

target_range = init_range + (target_speed * sim_time(:));  % time varying target range
target_range = reshape(target_range, n_samples, n_pulses); % restore shape n_samples x n_pulses
time_delay = 2*target_range/c;                             % time delay associated to target range

delta_time   = (fast_time * ones(1, n_pulses)) - time_delay;               % array time - time back and forth to target
rx_signal    = exp(1i*2*pi*(fc*delta_time + slope * ((delta_time).^2)/2)); % received signal

% generate mixed signal
mixed_signal = (tx_signal* ones(1, n_pulses)).* conj(rx_signal); % mix tranmitted and received signals to obtain beat signal

% RDM axes
rangeBinAxis = (0:n_samples-1).*c/(2*bandwidth);
dopplerBinSize = (1/chirp_duration)/n_pulses;
velocityBinAxis = (-n_pulses/2:n_pulses/2-1).*dopplerBinSize*lambda/2;

% 2D FFT to perform range and Doppler map
range_doppler_map = fftshift(fft2(mixed_signal), 2);
range_doppler_map = range_doppler_map/(max(max(abs(range_doppler_map))));

% Plot the RDM for the valid ranges of interest - targets ahead of you
figure(1);
imagesc(velocityBinAxis, rangeBinAxis, 20*log10(abs(range_doppler_map)));
colorbar;
ylabel("Range (m)");
xlabel("Velocity (m/s)");

