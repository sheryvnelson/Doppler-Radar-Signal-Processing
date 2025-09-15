%% Radar Doppler Demo - Compare Towards vs Away Motion
% - Plots the I/Q signals and FFTs together to compare
% - Estimates speed using Doppler principle

clc; clear; close all;

%% --- Setup Radar ---
addpath('..\\..\\RadarSystem'); % Path to RadarSystem library
resetRS;
szPort = findRSPort;
if isempty(szPort)
error('Radar device not found.');
end
oRS = RadarSystem(szPort);

% Radar parameters
NTS = 256; % Number of samples per frame
frame_time = 0.15; % Frame time in seconds
Nfft = 512; % FFT length

% Apply parameters to radar
oRS.oEPRadarS2GLP.parameters.number_of_samples = NTS;
oRS.oEPRadarS2GLP.parameters.frame_time_sec = frame_time;
oRS.oEPRadarS2GLP.apply_parameters;

% Radar constants
c = 3e8; % Speed of light
f_c = 24.005e9; % Radar center frequency (24.005 GHz)
lambda = c / f_c; % Wavelength
Fs = NTS / frame_time; % Sampling frequency

disp('Radar connected and configured.');
disp('----------------------------------------');

%% --- Step 1: Capture frame while moving TOWARDS radar ---
disp('> Please get ready to move TOWARDS the radar.');
disp('> When ready, press Enter. You will have 0.5 sec to start moving.');
pause; % Wait for user
disp('Capturing frame in...');
pause(0.5); % Give user time to move

raw1 = oRS.oEPRadarS2GLP.get_raw_data;
frame_towards = double(raw1.sample_data(:,1));
disp('Frame captured while moving towards.');

%% --- Step 2: Capture frame while moving AWAY from radar ---
disp('> Now get ready to move AWAY from the radar.');
disp('> When ready, press Enter. You will have 0.5 sec to start moving.');
pause;
disp('Capturing frame in...');
pause(0.5); % Give user time to move

raw2 = oRS.oEPRadarS2GLP.get_raw_data;
frame_away = double(raw2.sample_data(:,1));
disp('Frame captured while moving away.');

%% --- Step 3: Process and plot ---
% Apply Hamming window to reduce spectral leakage
win = hamming(NTS);

% Compute FFT (frequency domain)
fft_towards = fftshift(fft(frame_towards .* win, Nfft));
fft_away = fftshift(fft(frame_away .* win, Nfft));

% Frequency axis (Hz)
freq = linspace(-Fs/2, Fs/2, Nfft);

% DC suppression
dc_range = round(0.02 * Nfft); % 2% of spectrum
center = round(Nfft/2);
fft_towards(center - dc_range : center + dc_range) = 0;
fft_away(center - dc_range : center + dc_range) = 0;

%% --- Step 4: Plot I/Q signals ---
figure('Name','I/Q Comparison','NumberTitle','off');
plot(real(frame_towards), 'b','DisplayName','Towards'); hold on;
plot(real(frame_away), 'r','DisplayName','Away');
title('I Component: Towards (blue) vs Away (red)');
xlabel('Sample Index'); ylabel('Amplitude');
legend; grid on;


%% --- Step 5: Plot FFT magnitude spectra ---
figure('Name','FFT Comparison','NumberTitle','off');
plot(freq/1e3, abs(fft_towards), 'b','DisplayName','Towards'); hold on;
plot(freq/1e3, abs(fft_away), 'r','DisplayName','Away');
title('FFT Magnitude: Towards (blue) vs Away (red)');
xlabel('Frequency (kHz)'); ylabel('Magnitude');
legend; grid on;

%% --- Step 6: Estimate speed ---
% Find peak in magnitude spectrum
[~, idx_towards] = max(abs(fft_towards));
[~, idx_away] = max(abs(fft_away));

% Doppler frequency
f_d_towards = freq(idx_towards);
f_d_away = freq(idx_away);


% Doppler speed calculation: v = (f_d * lambda) / 2
speed_towards = (f_d_towards * lambda) / 2;
speed_away = (f_d_away * lambda) / 2;

disp('----------------------------------------');
disp(['Estimated speed while moving towards: ' num2str(speed_towards, '%.2f') ' m/s']);
disp(['Estimated speed while moving away: ' num2str(speed_away, '%.2f') ' m/s']);
disp('----------------------------------------');




disp('Demo finished. Compare plots and speed values.');




