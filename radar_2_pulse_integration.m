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
n_integration = 5;                                % number of pulse integration

% time vector
n_samples  = 2 * max_range / range_resolution;        % number of samples
oversample = 10;                                      % oversample for smooth representation
t = linspace(0, 2*max_range/c, oversample*n_samples); % time vector

% generate a rectangular pulse
tx_pulse = rectpuls(t - pulse_duration, pulse_duration);

% define target parameters
target_range = 10e3;                 % [m] target range
target_rcs = 1;                      % [m²] target radar cross section
target_delay = 2 * target_range / c; % [s] round trip delay

% simulate the received pulse with delay
rx_pulse = target_rcs * rectpuls(t - pulse_duration - target_delay, pulse_duration);
rx_pulse = ones(n_integration,1) * rx_pulse; % duplicate received signal n_integration times
% add white Gaussian noise to the received signal
SNR = 0;                    % [dB]
noise_power = 10^(-SNR/10); % noise power
rng(2025);                  % fix random seed for reproductibility
rx_pulse_noisy = rx_pulse + sqrt(noise_power) * randn(size(rx_pulse));

% matched filter (correlation with transmitted pulse)
for i_integration = 1:n_integration
  mf_output(i_integration,:) = xcorr(rx_pulse_noisy(i_integration,:) , tx_pulse);
endfor
% perform non coherent integration
mf_output_integration = sum(abs(mf_output), 1);

% generate range axis
range_axis = linspace(-max_range, max_range, length(mf_output));

% CFAR implementation
n_training_cells = 10; % number of training cells
n_guard_cells = 5;     % number of guard cells

cfar_windows=(1/n_training_cells)*[ones(1,n_training_cells/2),zeros(1,n_guard_cells),ones(1,n_training_cells/2)];

mf_output_dwn = downsample(mf_output_integration, oversample);
range_axis_dwn = downsample(range_axis, oversample);
cfarThreshold = conv(mf_output_dwn, cfar_windows, 'same');
% TODO check offset formula - specially with pulse intergration
Pfa = 1e-6;
offset = log(n_training_cells)*(Pfa^(-1/(n_training_cells))-1);

% plot transmitted and received pulses
figure(1);
subplot(3,1,1);
plot(t, rx_pulse_noisy(1,:));
title('Single Received Pulse with Noise');
xlabel('Time [s]');
ylabel('Amplitude');
grid('on');

% plot matched filter output
subplot(3,1,2);
plot(range_axis, abs(mf_output(1,:)));
title('Single Matched Filter Output');
xlabel('Range [m]');
ylabel('Amplitude');
grid('on');
xlim([0 max_range]);

subplot(3,1,3);
plot(range_axis, abs(mf_output(1,:)));
title('Single Matched Filter Output - Zoom');
xlabel('Range [m]');
ylabel('Amplitude');
grid('on');
xlim([(target_range - 2e3) (target_range + 2e3)]);

figure(3);
plot(range_axis_dwn, mf_output_dwn, range_axis_dwn, cfarThreshold + offset,'r--');
xlim([0 max_range]);
title('CA-CFAR threshold with pulse integration');
xlabel('Range [m]');
ylabel('Amplitude');
legend({'match filter output','CA CFAR threshold'},'Location','northeast');
grid('on');
