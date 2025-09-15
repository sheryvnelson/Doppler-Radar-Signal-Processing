clc; clear; close all;

% Radar Setup
addpath('..\\..\\RadarSystem');
resetRS;
szPort = findRSPort;
if isempty(szPort), error('Radar not found'); end
oRS = RadarSystem(szPort);

% Radar Parameters
NTS = 256; % Samples per frame
frame_time = 0.15; % Frame time
Nfft = 512;
numFrames = 10;

oRS.oEPRadarS2GLP.parameters.number_of_samples = NTS;
oRS.oEPRadarS2GLP.parameters.frame_time_sec = frame_time;
oRS.oEPRadarS2GLP.apply_parameters;

% Radar Constants
Fs = 1 / (frame_time / NTS);
c = 3e8;
f_c = 24.005e9;
lambda = c / f_c;

% Preallocate
storedFrames = complex(zeros(NTS, numFrames));
estimatedSpeeds = zeros(1, numFrames);

disp('Collecting 10 radar frames..');

%frames
for i = 1:numFrames
pause(frame_time);
raw = oRS.oEPRadarS2GLP.get_raw_data;
iq = double(raw.sample_data(:,1));
storedFrames(:, i) = iq;
end

% Process and plot each frame
for i = 1:numFrames
iq = storedFrames(:, i);
I = real(iq);
Q = imag(iq);

% Windowed I/Q
window = hamming(NTS);
iqWin = iq .* window;

% FFT
fftData = fftshift(fft(iqWin, Nfft));
mag = 20*log10(abs(fftData) + eps);
freqAxis = linspace(-Fs/2, Fs/2, Nfft)/1e3; % kHz

% DC suppression
dcCenter = round(Nfft/2);
dcRange = round(0.02 * Nfft);
mag(dcCenter-dcRange:dcCenter+dcRange) = -100;

% Estimate speed 
[~, peakIdx] = max(mag);
f_d = freqAxis(peakIdx) * 1e3; % Doppler frequency in Hz 
speed = (f_d * lambda) / 2; % Speed in m/s 
estimatedSpeeds(i) = speed;

% Plot figure
figure('Name', ['Frame ' num2str(i)], 'NumberTitle', 'off', 'Position', [100 100 1200 700]);

subplot(2,2,1);
plot(I, 'b'); hold on;
plot(Q, 'r');
title(['Raw I/Q - Frame ' num2str(i)]);
xlabel('Sample'); ylabel('Amplitude'); legend('I','Q');
grid on;

subplot(2,2,2);
plot(real(iqWin), 'b'); hold on;
plot(imag(iqWin), 'r');
title('Windowed I/Q');
xlabel('Sample'); ylabel('Amplitude'); legend('I','Q');
grid on;

subplot(2,2,3);
plot(freqAxis, mag, 'b');
title('FFT Magnitude Spectrum');
xlabel('Frequency (kHz)'); ylabel('Magnitude (dB)');
grid on;

subplot(2,2,4);
plot(1:i, estimatedSpeeds(1:i), '-o', 'Color', 'm', 'LineWidth', 2);
title('Estimated Speed Over Frames');
xlabel('Frame'); ylabel('Speed (m/s)');
grid on;
xlim([1 numFrames]);
ylim([-max(abs(estimatedSpeeds)) - 0.5, max(abs(estimatedSpeeds)) + 0.5]);
fprintf('Frame %d - Estimated Speed: %.2f m/s\n', i, speed);
end

disp('✅ All 10 frames processed and plotted with signed speed.');