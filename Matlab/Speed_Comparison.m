%% ================== Live Doppler: 5 Frames with Delay ==================

clc; clear; close all;

%% --- Setup Radar  demo ---
addpath('..\\..\\RadarSystem'); % Path to RadarSystem library
resetRS;
szPort = findRSPort;
if isempty(szPort)
    error('Radar device not found.');
end
oRS = RadarSystem(szPort);

% Radar parameters
NTS        = 256;     % samples per frame
frame_time = 0.15;    % seconds per frame 
Nfft       = 512;     % FFT length 

% Apply parameters to radar
oRS.oEPRadarS2GLP.parameters.number_of_samples = NTS;
oRS.oEPRadarS2GLP.parameters.frame_time_sec    = frame_time;
oRS.oEPRadarS2GLP.apply_parameters;

% Radar constants
c     = 3e8;
f_c   = 24.005e9;     
lambda= c / f_c;
Fs    = NTS / frame_time; 

disp('Radar connected and configured.');
disp('----------------------------------------');

numFramesToProcess = 5;    
frameDelay         = 0.75; 
dcFrac             = 0.02; 
labelThresh        = 0.05; 

%% --- Precompute ---
assert(mod(Nfft,2)==0, 'Nfft must be even.');
dcRange    = max(1, round(dcFrac * Nfft));
labels     = strings(1, numFramesToProcess);
estSpeeds  = zeros(1, numFramesToProcess);

win        = hamming(NTS);
freqAxisHz = linspace(-Fs/2, Fs/2, Nfft);   
freqAxiskHz= freqAxisHz / 1e3;              
[~, dcCenter] = min(abs(freqAxisHz));       

%% --- Figure (reused) ---
fig = figure('Name','Live Doppler (5 frames)','NumberTitle','off','Position',[100 100 1200 700]);

disp('Starting capture... move however you like (towards/away).');

for i = 1:numFramesToProcess
   
    % get_raw_data returns struct with .sample_data
    raw = oRS.oEPRadarS2GLP.get_raw_data;
    sd  = double(raw.sample_data); 
    % Case A: complex in one column  
    if ~isreal(sd(:,1))
        iq = sd(:,1);
    % Case B: two real columns [I, Q]
    elseif size(sd,2) >= 2
        iq = complex(sd(:,1), sd(:,2));
    % Case C: only real column available → treat as real (imag=0)
    else
        iq = complex(sd(:,1), 0);
    end

    if numel(iq) ~= NTS
        warning('Received %d samples, expected %d. Truncating/padding.', numel(iq), NTS);
        iq = iq(:);
        if numel(iq) >= NTS
            iq = iq(1:NTS);
        else
            iq = [iq; zeros(NTS-numel(iq),1)];
        end
    end

    I = real(iq); Q = imag(iq);

    %% --------- PROCESS ---------
    iqWin   = iq .* win;
    fftData = fftshift(fft(iqWin, Nfft));
    magdB   = 20*log10(abs(fftData) + eps);

    % DC suppression 
    lo = max(1, dcCenter - dcRange);
    hi = min(Nfft, dcCenter + dcRange);
    magdB(lo:hi) = -100;

    % Peak search 
    search = magdB; search(~isfinite(search)) = -Inf;
    [~, peakIdx] = max(search);

    f_d   = freqAxisHz(peakIdx);   % Doppler frequency (Hz)
    speed = (f_d * lambda) / 2;    % m/s (CW Doppler)
    estSpeeds(i) = speed;

    % Motion label
    if abs(speed) < labelThresh
        labels(i) = "Stationary";
    elseif speed > 0
        labels(i) = "Approaching";
    else
        labels(i) = "Receding";
    end

    %% --------- PLOT (reuse figure) ---------
    clf(fig);

    % Raw I/Q
    subplot(2,2,1);
    plot(I, 'b'); hold on; plot(Q, 'r');
    title(sprintf('Raw I/Q - Frame %d', i));
    xlabel('Sample'); ylabel('Amplitude'); legend('I','Q'); grid on;

    % Windowed I/Q
    subplot(2,2,2);
    plot(real(iqWin), 'b'); hold on; plot(imag(iqWin), 'r');
    title('Windowed I/Q'); xlabel('Sample'); ylabel('Amplitude'); legend('I','Q'); grid on;

    % FFT Magnitude with annotation
    subplot(2,2,3);
    plot(freqAxiskHz, magdB, 'b'); hold on;
    plot(freqAxiskHz(peakIdx), magdB(peakIdx), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    title('FFT Magnitude Spectrum');
    xlabel('Frequency (kHz)'); ylabel('Magnitude (dB)'); grid on;
    text(freqAxiskHz(peakIdx), magdB(peakIdx)+3, ...
        sprintf('f_d = %.1f Hz', f_d), ...
        'Color','k','FontSize',10,'FontWeight','bold');

    % Speed history
    subplot(2,2,4);
    plot(1:i, estSpeeds(1:i), '-o', 'LineWidth', 2); hold on;
    plot(i, estSpeeds(i), 'ko', 'MarkerFaceColor','y', 'MarkerSize', 8);
    title('Estimated Speed Over Frames');
    xlabel('Frame'); ylabel('Speed (m/s)'); grid on;
    xlim([1 numFramesToProcess]);
    maxRange = max(1, max(abs(estSpeeds(1:i))) + 0.5);
    ylim([-maxRange, maxRange]);
    text(i, estSpeeds(i), sprintf('%.2f m/s (%s)', speed, labels(i)), ...
        'VerticalAlignment','bottom','HorizontalAlignment','right', ...
        'FontSize',10,'FontWeight','bold','BackgroundColor','w');

    drawnow limitrate;

    %% --------- LOG ---------
    fprintf('Frame %d:\n', i);
    fprintf(' → Doppler Frequency (f_d): %.2f Hz\n', f_d);
    fprintf(' → Estimated Speed: %.2f m/s (%s)\n', speed, labels(i));
    fprintf('----------------------------------------\n');

    pause(frameDelay);
end

disp('✅ 5 frames captured and processed with fresh radar data.');
