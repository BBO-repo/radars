clc;
clear all;
close all;

if ( exist('OCTAVE_VERSION', 'builtin') ~= 0 )
  % load signal package
  run .octaverc
endif

%% radar design parameters
% Pmin = Pt[G²λ²σ]/[(4π)³r⁴] % [W] minimum detectable power
% could not find order of magnitude of min detectable power so set max range instead
max_range = 30e3;      % [m] maximum range
pulse_duration = 1e-6; % [s] pulse time width

%% determine other radar parameters according to design parameters
c = 3e8;                                          % [m/s] speed of light
pulse_repetition_frequency = c / (2 * max_range); % [Hz] pulse repetition frequency (Hz)
range_resolution = c * pulse_duration / 2;        % [m] range resolution
fc = 10e9;                                        % [Hz] carrier frequency

% time vector
n_samples  = 2 * max_range / range_resolution;        % number of samples
oversample = 10;                                      % oversample for smooth representation
t = linspace(0, 2*max_range/c, oversample*n_samples); % time vector

% generate a rectangular pulse
tx_pulse = rectpuls(t - pulse_duration, pulse_duration);

% transmitted signal spectrum
NFFT = 2^nextpow2(length(t));                      % number of FFT bins
fs = 1/(t(2)-t(1));                                % sample frequency
f = (fs/NFFT) * (0:NFFT-1) - 0.5*fs;               % frequency vector
tx_pulse_spectrum = fftshift(fft(tx_pulse, NFFT)); % compute the FFT

% define target parameters
target_range = 8e3;                  % [m] target range
target_rcs = 1;                      % [m²] target radar cross section
target_delay = 2 * target_range / c; % [s] round trip delay

% simulate the received pulse with delay
rx_pulse = target_rcs * rectpuls(t - pulse_duration - target_delay, pulse_duration);

% add white Gaussian noise to the received signal
SNR = 10;                   % [dB]
noise_power = 10^(-SNR/10); % noise power
rng(2025);                  % fix random seed for reproductibility
rx_pulse_noisy = rx_pulse + sqrt(noise_power) * randn(size(rx_pulse));

% matched filter (correlation with transmitted pulse)
mf_output = xcorr(rx_pulse_noisy, tx_pulse);

% generate range axis
range_axis = linspace(-max_range, max_range, length(mf_output));

% CFAR implementation
n_training_cells = 10; % number of training cells
n_guard_cells = 5;     % number of guard cells

cfar_windows=(1/n_training_cells)*[ones(1,n_training_cells/2),zeros(1,n_guard_cells),ones(1,n_training_cells/2)];

mf_output_dwn = downsample(mf_output, oversample);
range_axis_dwn = downsample(range_axis, oversample);
cfarThreshold = conv(mf_output_dwn, cfar_windows, 'same');
% TODO check offset formula
Pfa = 1e-5;                                                     % false alarm rate
offset = log(n_training_cells)*(Pfa^(-1/(n_training_cells))-1); % cfar offset

% detected targets
idx = find(mf_output_dwn > (cfarThreshold + offset) );
display(sprintf('target detected at %d[m] \n',range_axis_dwn(idx)));

% plot transmitted and received pulses
figure(1);
subplot(3,1,1);
plot(t, tx_pulse);
title('Transmitted Pulse');
xlabel('Time [s]');
ylabel('Amplitude');
grid('on');

subplot(3,1,2);
stem(t, tx_pulse,'LineStyle','none');
title('Transmitted Pulse - Zoom');
xlabel('Time [s]');
ylabel('Amplitude');
grid('on');
xlim([0 2*pulse_duration]);

subplot(3,1,3);
plot(f,abs(tx_pulse_spectrum));
title('Transmitted Pulse Spectrum');
xlabel('Frequency [Hz]');
ylabel('Amplitude');
grid('on');

figure(2);
subplot(3,1,1);
plot(t, rx_pulse_noisy);
title('Received Pulse with Noise');
xlabel('Time [s]');
ylabel('Amplitude');
grid('on');

% plot matched filter output
subplot(3,1,2);
plot(range_axis, abs(mf_output));
title('Matched Filter Output');
xlabel('Range [m]');
ylabel('Amplitude');
grid('on');
xlim([0 max_range]);

subplot(3,1,3);
plot(range_axis, abs(mf_output));
title('Matched Filter Output - Zoom');
xlabel('Range [m]');
ylabel('Amplitude');
grid('on');
xlim([(target_range - 2e3) (target_range + 2e3)]);

figure(3);
plot(range_axis_dwn, mf_output_dwn, range_axis_dwn, cfarThreshold + offset,'r--');
xlim([0 max_range]);
title('CA-CFAR threshold');
xlabel('Range [m]');
ylabel('Amplitude');
legend({'match filter output','CA CFAR threshold'},'Location','northeast');
grid('on');
