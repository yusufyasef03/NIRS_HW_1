%% NIRS Homework 1: MATLAB Version (Two Channels)
clear all; close all; clc;

% 1. LOAD DATA
fprintf('Loading data files...\n');
% Ensure the file is in your current MATLAB folder
load('FN_031_V2_Postdose2_Nback.mat'); 

% 2. CONVERT INTENSITY TO OPTICAL DENSITY (OD)
% Calculate change in OD relative to the mean
baseline = mean(d);
dc = -log(d ./ baseline); 

[nTime, nChan] = size(d);
half = nChan / 2; % Usually 690nm and 830nm

% 3. PLOT RAW DATA & OD (Only first 2 channels)
figure(1);
subplot(2,1,1);
plot(t, d(:, 1:2), 'LineWidth', 1); 
title('Raw Intensity (Channels 1 and 2)');
xlabel('Time (s)'); ylabel('Intensity');
legend('Ch1', 'Ch2');

subplot(2,1,2);
plot(t, dc(:, 1:2), 'LineWidth', 1); 
title('Change in Optical Density (\DeltaOD)');
xlabel('Time (s)'); ylabel('\DeltaOD');
legend('Ch1', 'Ch2');

% 4. BEER-LAMBERT LAW (Calculate HbO and HbR)
% Standard extinction coefficients
E = [0.0955, 0.4323;   
     0.2526, 0.1798];
invE = inv(E);

% Path length: SDS (3cm) * DPF (6)
dist = 3 * 6; 

HbO = zeros(nTime, half);
HbR = zeros(nTime, half);

% Calculate concentration for every channel
for i = 1:half
    temp_OD = [dc(:, i), dc(:, i+half)]'; 
    conc = invE * (temp_OD / dist);
    
    HbO(:, i) = conc(1, :)'; 
    HbR(:, i) = conc(2, :)'; 
end

% 5. FILTERING (Basic MATLAB Moving Average)
% We use a window of 15 samples to smooth the noise
windowSize = 15; 
b_coeff = (1/windowSize) * ones(1, windowSize);

HbO_filt = filter(b_coeff, 1, HbO);
HbR_filt = filter(b_coeff, 1, HbR);

% 6. FOURIER TRANSFORM (0.5Hz to 2.5Hz)
fs = 1 / (t(2) - t(1)); 
L = length(d(:,1));
Y = fft(d(:,1));
P2 = abs(Y/L);
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(floor(L/2)))/L;

figure(2);
plot(f, P1, 'Color', [0.5 0 0.5], 'LineWidth', 1.5);
title('Frequency Spectrum (0.5 - 2.5 Hz)');
xlabel('Frequency (Hz)'); ylabel('|P1(f)|');
xlim([0.5 2.5]); 
grid on;

% 7. FINAL RESULTS PLOT (Channel 1 and Channel 2)
figure(3);

% Subplot for Channel 1
subplot(2,1,1);
plot(t, HbO_filt(:, 1), 'r', 'LineWidth', 1.5); hold on;
plot(t, HbR_filt(:, 1), 'b', 'LineWidth', 1.5);
title('Concentration Changes - Channel 1');
ylabel('\Delta Conc (\muM)');
legend('HbO', 'HbR');
grid on;

% Subplot for Channel 2
subplot(2,1,2);
plot(t, HbO_filt(:, 2), 'r', 'LineWidth', 1.5); hold on;
plot(t, HbR_filt(:, 2), 'b', 'LineWidth', 1.5);
title('Concentration Changes - Channel 2');
xlabel('Time (s)'); ylabel('\Delta Conc (\muM)');
legend('HbO', 'HbR');
grid on;
