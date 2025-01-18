clc;
clear all;
close all;

if ( exist('OCTAVE_VERSION', 'builtin') ~= 0 )
  % load signal package
  run .octaverc
endif

%% radar design parameters
pulse_duration = 10e-6;           % [s] pulse time width
bandwidth = 10e6;                 % [Hz] LFM sweep bandwidth
slope = bandwidth/pulse_duration; % [Hz/s] slope of the linear LFM
fc = 10e9;                        % [Hz] carrier frequency

%% determine other radar parameters according to design parameters
n_samples  = 4096;   % number of samples
fs = 100e6;          % sample frequency
c  = 3e8;            % [m/s] speed of light

% time and frequency vector
t  = (0:n_samples-1)/fs;
range_axis = c/2 .* t;
f = (-n_samples/2:n_samples/2-1)*fs/n_samples;

% generate a rectangular LFM pulse
tx_pulse      = rectpuls(t-pulse_duration/2,pulse_duration).*exp(1i*pi*slope*(t-pulse_duration/2).^2);
tx_pulse_hamm = rectpuls(t-pulse_duration/2,pulse_duration).*exp(1i*pi*slope*(t-pulse_duration/2).^2).*[hamming(1001)',zeros(1,(n_samples-1001))];

% define target parameters
target_range = 4e3;
target_delay = 2*target_range/c;

% simulate the received pulse with delay
rx_pulse   = rectpuls(t-pulse_duration/2-target_delay,pulse_duration).*exp(1i*pi*slope*(t-pulse_duration/2-target_delay).^2).*exp(-1i*2*pi*fc*target_delay);

% add white Gaussian noise to the received signal
SNR = 20;                   % [dB]
noise_power = 10^(-SNR/10); % noise power
rng(2025);                  % fix random seed for reproductibility
rx_pulse = rx_pulse + sqrt(noise_power) * randn(size(rx_pulse));

%FFT
rx_pulse_spectrum = fft(rx_pulse, n_samples);
tx_pulse_spectrum = fft(tx_pulse, n_samples);

% matched filter
tx_pulse_hamm_spectrum = fft(tx_pulse_hamm,n_samples);
tx_pulse_spectrum = fft(tx_pulse,n_samples);
mf_output = ifft(rx_pulse_spectrum.*conj(tx_pulse_hamm_spectrum));
##mf_output = ifft(rx_pulse_noisy.*conj(tx_pulse_spectrum));


% plot transmitted and received pulses
figure(1);
subplot(2,1,1);
plot(range_axis,real(tx_pulse));
grid('on');
xlabel('Range [m]');
ylabel('Amplitude');
title('Transmitted Pulse (Re)');

subplot(2,1,2);
plot(f,abs(fftshift(tx_pulse_spectrum)));
grid('on');
xlabel('Frequency [Hz]');
ylabel('Amplitude');
title('Transmitted Pulse Spectrum');

figure(2);
subplot(2,1,1);
plot(range_axis,real(rx_pulse));
grid('on');
xlabel('Range [m]');
ylabel('Amplitude');
title('Received Pulse with Noise (Re)');

subplot(2,1,2);
plot(f,abs(fftshift(rx_pulse_spectrum)));
grid('on');
xlabel('Frequency [Hz]');
ylabel('Amplitude');
title('Spectrum of signal');

% Pulse Compression
figure(3);
plot(range_axis,20*log10(abs(mf_output)/max(abs(mf_output))));
grid('on');
xlabel('Range [m]');
ylabel('Amplitude [dB]');
title('Result of pulse compression');
axis('tight');

