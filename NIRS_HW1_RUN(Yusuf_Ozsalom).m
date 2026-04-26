clear all; close all; clc;

% --- TASK 1 & 2: Run MBLL Conversion ---
% Set basic params
SDS = 3; 
tissueType = 'adult_head'; 
plotChannels = [1, 2]; 

% Process Subject 1
fprintf('Running Subject 1...\n');
[dHbR_1, dHbO_1, fig1] = CalcNIRS('FN_031_V2_Postdose2_Nback.mat', SDS, ...
    tissueType, plotChannels, 'ExtinctionCoefficientsData.csv', ...
    'DPFperTissue.txt', 'RelativeDPFCoefficients.csv');

% Process Subject 2
fprintf('Running Subject 2...\n');
[dHbR_2, dHbO_2, fig2] = CalcNIRS('FN_032_V1_Postdose1_Nback.mat', SDS, ...
    tissueType, plotChannels, 'ExtinctionCoefficientsData.csv', ...
    ... % (keeping file inputs same as above)
    'DPFperTissue.txt', 'RelativeDPFCoefficients.csv');

% --- TASK 3: FFT & SNR (Check Signal Quality) ---
% Pull sampling rate and signal length
data = load('FN_031_V2_Postdose2_Nback.mat');
Fs = 1 / mean(diff(data.t)); 
L = length(dHbO_1(:,1)); 

% Math to move from Time -> Frequency domain
Y = fft(dHbO_1(:, 1));
P2 = abs(Y/L);
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:floor(L/2))/L; 

% Plot it (Zooming to 0-5Hz to see body signals)
figure;
plot(f, P1, 'r', 'LineWidth', 1.5);
title('Task 3: FFT (Ch1, Sub1)');
xlabel('Freq (Hz)'); ylabel('Amp');
xlim([0 5]); 

% Calculate SNR: Heartbeat peak vs. background noise floor
noise_idx = find(f > 2.5);
noise = mean(P1(noise_idx));

hb_idx = find(f > 0.8 & f < 2.0); % Range where pulse usually sits
signal = max(P1(hb_idx));

SNR = signal / noise;

% Print stats
fprintf('\n--- TASK 3 --- \n');
fprintf('Pulse Peak: %f\n', signal);
fprintf('Noise Floor: %f\n', noise);
fprintf('SNR: %f\n', SNR);
