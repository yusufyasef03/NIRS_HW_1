%% NIRS Homework 1: MATLAB Version
clear all; close all; clc;

% 1. LOAD DATA
% We load the subject file. Note: pkg load is NOT used here for MATLAB.
fprintf('Loading data files...\n');
if exist('FN_031_V2_Postdose2_Nback.mat', 'file')
    load('FN_031_V2_Postdose2_Nback.mat'); 
else
    error('Data file not found. Please check the file name.');
end

% 2. CONVERT INTENSITY TO OPTICAL DENSITY (OD)
% We calculate the change in OD (dc) relative to the mean.
baseline = mean(d);
dc = -log(d ./ baseline); 

[nTime, nChan] = size(d);
half = nChan / 2; % 690nm and 830nm channels

% 3. PLOT RAW DATA & OD (Task: Only first 2 channels)
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
% Standard extinction coefficients for 690nm and 830nm.
E = [0.0955, 0.4323;   
     0.2526, 0.1798];
invE = inv(E);

% Path length: SDS (3cm) * DPF (6).
dist = 3 * 6; 

HbO = zeros(nTime, half);
HbR = zeros(nTime, half);

% Formula: [HbO; HbR] = E^-1 * (Delta OD / Distance)
for i = 1:half
    temp_OD = [dc(:, i), dc(:, i+half)]'; 
    conc = invE * (temp_OD / dist);
    
    HbO(:, i) = conc(1, :)'; 
    HbR(:, i) = conc(2, :)'; 
end

% 5. FILTERING (MATLAB Basic Version)
% We use a Moving Average filter because 'butter' requires a toolbox.
% This removes high-frequency noise like the heart rate.
windowSize = 15; 
b_coeff = (1/windowSize) * ones(1, windowSize);

HbO_filt = filter(b_coeff, 1, HbO);
HbR_filt = filter(b_coeff, 1, HbR);

% 6. FOURIER TRANSFORM (Frequency Analysis)
fs = 1 / (t(2) - t(1)); % Sample Rate
L = length(d(:,1));
Y = fft(d(:,1));
P2 = abs(Y/L);
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(floor(L/2)))/L;

figure(2);
plot(f, P1, 'Color', [0.5 0 0.5], 'LineWidth', 1.5);
title('Frequency Spectrum of Raw Data');
xlabel('Frequency (Hz)'); ylabel('|P1(f)|');
% Focus on 0.5Hz to 2.5Hz to see the heart rate peak
xlim([0.5 2.5]); 
grid on;

% 7. FINAL RESULTS PLOT
figure(3);
plot(t, HbO_filt(:, 1), 'r', 'LineWidth', 1.5); hold on;
plot(t, HbR_filt(:, 1), 'b', 'LineWidth', 1.5);
title('Concentration Changes (Moving Average) - Channel 1');
xlabel('Time (s)'); ylabel('\Delta Conc (\muM)');
legend('HbO (Oxy)', 'HbR (Deoxy)');
grid on;
