clear all; close all; clc;

% Basic setup based on the lab manual
SDS = 3; 
plotChannels = [1, 2]; % We only need to check the first two for the HW

% Hardcoded filenames so we don't have to type them every time
file1 = 'FN_031_V2_Postdose2_Nback.mat';
file2 = 'FN_032_V1_Postdose1_Nback.mat';
ext   = 'ExtinctionCoefficientsData.csv';
dpf   = 'DPFperTissue.txt';
rel   = 'RelativeDPFCoefficients.csv';

% --- Task 1 & 2: Get the Hemoglobin Data ---
fprintf('Starting Subject 1 Analysis...\n');
[dHbR_1, dHbO_1, fig1] = CalcNIRS(file1, SDS, 'adult_head', plotChannels, ext, dpf, rel);

fprintf('Starting Subject 2 Analysis...\n');
[dHbR_2, dHbO_2, fig2] = CalcNIRS(file2, SDS, 'adult_head', plotChannels, ext, dpf, rel);

% --- Task 3: Signal Quality (FFT) ---
data = load(file1);
Fs = 1 / mean(diff(data.t)); % Calculate how fast data was recorded
L = length(dHbO_1(:,1));

% Run the math to switch from time to frequency
Y = fft(dHbO_1(:,1));
P1 = abs(Y/L); 
P1 = P1(1:floor(L/2)+1); 
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:floor(L/2))/L; % Frequency axis in Hz

% Calculate SNR (Heartbeat vs Background noise)
noise = mean(P1(f > 2.5)); % Pure noise area
signal = max(P1(f > 0.8 & f < 2.0)); % Where the heart pulse lives
SNR = signal / noise;

fprintf('\n--- TASK 3 STATS ---\n');
fprintf('Subject 1 SNR: %.4f\n', SNR);

% Plot the FFT to see the heartbeat peak
figure; plot(f, P1, 'r'); xlim([0 5]);
title('Frequency Check (Look for Heartbeat)');
xlabel('Hz'); ylabel('Amplitude');
