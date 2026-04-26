%% NIRS Homework 1: Data Processing
clear all; close all; clc;

% 1. LOAD EXPERIMENTAL DATA
% Loading the first subject's .mat file
fprintf('Loading subject data...\n');
load('FN_031_V2_Postdose2_Nback.mat'); 

% 2. LOAD AUXILIARY FILES
% Using 'importdata' to handle mixed text/numeric formats
fprintf('Loading coefficients...\n');
ext_coeffs = importdata('ExtinctionCoefficientsData.csv');
rel_dpf = importdata('RelativeDPFCoefficients.csv');

% Handle the problematic TXT file
% We use 'fileread' to see it or 'importdata'
try
    dpf_info = importdata('DPFperTissue.txt');
    disp('DPF File Loaded successfully.');
catch
    fprintf('Warning: Manual check of DPFperTissue.txt may be needed.\n');
end

% 3. CONVERSION: INTENSITY TO OPTICAL DENSITY (OD)
% 'd' is your raw intensity matrix. 
% We normalize by the mean of each channel to get change in OD (dc)
baseline = mean(d);
dc = -log(d ./ baseline); 

% 4. DATA DIMENSIONS
[nTimePoints, nChannels] = size(d);
halfChans = nChannels / 2; % Usually NIRS has 2 wavelengths (e.g. 690nm & 830nm)

% 5. VISUALIZATION (Task 1: Raw Data Check)
figure(1);
subplot(2,1,1);
plot(t, d(:, 1:5)); % Plot first 5 channels of raw intensity
title('Raw Intensity (First 5 Channels)');
xlabel('Time (s)'); ylabel('Intensity');

subplot(2,1,2);
plot(t, dc(:, 1:5)); % Plot first 5 channels of Delta OD
title('Change in Optical Density (\DeltaOD)');
xlabel('Time (s)'); ylabel('\DeltaOD');

% 6. STIMULUS CHECK
% The 's' variable shows when the N-Back task occurred
figure(2);
plot(t, s);
title('Stimulus Design (N-Back Task Triggers)');
xlabel('Time (s)'); ylabel('On/Off');
ylim([-0.1 1.1]);

% 7. FORMATTING FOR NEXT STEP
fprintf('Data processed. Ready for Beer-Lambert Law application.\n');


% Find the indices where stimulus 's' is not zero
stim_indices = find(s > 0);
fprintf('Triggers found at time points: ');
disp(t(stim_indices)');

% Improved Stimulus Plot
figure(2);
stem(t, s, 'Marker', 'none', 'LineWidth', 1.5); 
title('Stimulus Design (N-Back Task Triggers)');
xlabel('Time (s)'); ylabel('Trigger (On/Off)');
grid on;
axis([0 max(t) -0.1 1.2]); % Force the axes to stay visible
