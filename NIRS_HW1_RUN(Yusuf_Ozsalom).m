clear all; close all; clc;

% Octave prep (signal toolbox for butter/filtfilt)
pkg load signal;

% 1. DATA IN
fprintf('Loading raw file...\n');
load('FN_031_V2_Postdose2_Nback.mat'); 

% 2. OPTICAL DENSITY (OD)
% Convert light intensity to absorption relative to average baseline
dc = -log(d ./ mean(d)); 
[nTime, nChan] = size(dc);
half = nChan / 2; % splitting 690nm and 830nm

% 3. MATH MAGIC (Beer-Lambert)
% Extinction coefficients (the "decoder")
E = [0.0955, 0.4323;   
     0.2526, 0.1798];
invE = inv(E);

% Path length: 3cm gap * DPF factor for scattering
L = 3; DPF = 6; 
dist = L * DPF;

HbO = zeros(nTime, half);
HbR = zeros(nTime, half);

% Crunch the pair-wise concentration
for i = 1:half
    % Formula: $[HbO; HbR] = E^{-1} \cdot (\Delta OD / Distance)$
    temp_OD = [dc(:, i), dc(:, i+half)]'; 
    conc = invE * (temp_OD / dist);
    
    HbO(:, i) = conc(1, :)';
    HbR(:, i) = conc(2, :)';
end

% 4. CLEAN THE JITTER
% Brain signals are slow (0.1Hz). Toss the heart rate/noise.
fs = 1 / (t(2) - t(1));  % Sample rate
fc = 0.1;                % LPF cutoff
[b, a] = butter(3, fc / (fs/2)); 

% Zero-phase filtering (no time-shift)
HbO_filt = filtfilt(b, a, HbO);
HbR_filt = filtfilt(b, a, HbR);

% 5. PLOTS
figure(1);
% Task timing (When the patient was actually thinking)
subplot(2,1,1);
stem(t, s, 'Marker', 'none', 'LineWidth', 1.5); 
title('N-Back Task Triggers');
grid on;

% The brain response (Channel 1)
subplot(2,1,2);
plot(t, HbO_filt(:, 1), 'r', 'LineWidth', 1.5, 'DisplayName', 'Oxy'); hold on;
plot(t, HbR_filt(:, 1), 'b', 'LineWidth', 1.5, 'DisplayName', 'Deoxy');
legend('show');
title('Filtered Signal - Chan 1');
xlabel('Time (s)'); 
ylabel('\Delta Conc (\muM)');
grid on;
axis tight;
