clear all; close all; clc;

% --- TASK 1 & 2: RUN THE FUNCTION ---
% Set the required inputs
SDS = 3; % 3cm assumed by PDF
tissueType = 'adult_head'; 
plotChannels = [1, 2]; % Plot first two channels as requested

% Run for Subject 1
fprintf('Running Subject 1...\n');
[dHbR_1, dHbO_1, fig1] = CalcNIRS('FN_031_V2_Postdose2_Nback.mat', SDS, tissueType, plotChannels, 'ExtinctionCoefficientsData.csv', 'DPFperTissue.txt', 'RelativeDPFCoefficients.csv');

% Run for Subject 2 (Required by PDF)
fprintf('Running Subject 2...\n');
[dHbR_2, dHbO_2, fig2] = CalcNIRS('FN_032_V1_Postdose1_Nback.mat', SDS, tissueType, plotChannels, 'ExtinctionCoefficientsData.csv', 'DPFperTissue.txt', 'RelativeDPFCoefficients.csv');

% --- TASK 3: FOURIER TRANSFORM & SNR ON FIRST FILE, FIRST CHANNEL ---
% Get time vector to calculate sampling frequency
data = load('FN_031_V2_Postdose2_Nback.mat');
Fs = 1 / mean(diff(data.t)); % Sampling frequency
L = length(dHbO_1(:,1));     % Length of signal

% Compute FFT
Y = fft(dHbO_1(:, 1));
P2 = abs(Y/L);
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:floor(L/2))/L; % Frequency vector

% Plot FFT
figure;
plot(f, P1, 'LineWidth', 1.5);
title('Task 3: Fourier Transform (Channel 1, Subject 1)');
xlabel('Frequency (Hz)'); ylabel('|P1(f)|');
xlim([0 5]); % Zoom in on 0 to 5 Hz to see the heartbeat clearly

% Calculate SNR
noise_idx = find(f > 2.5);
noise = mean(P1(noise_idx));

% Heartbeat is usually around 1 Hz to 1.5 Hz
hb_idx = find(f > 0.8 & f < 2.0); 
signal = max(P1(hb_idx));

SNR = signal / noise;
fprintf('\n--- TASK 3 RESULTS ---\n');
fprintf('Signal (Heartbeat peak): %f\n', signal);
fprintf('Noise (>2.5Hz average): %f\n', noise);
fprintf('SNR (Signal/Noise): %f\n', SNR);
