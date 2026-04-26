function [dHbR, dHbO, fig] = CalcNIRS(dataFile, SDS, tissueType, plotChannelIdx, extinctionCoefficientsFile, DPFperTissueFile, relDPFfile)
    
    % 1. Safety check: make sure the file actually exists before trying to open it
    if ~exist(dataFile, 'file')
        error('Data file does not exist.');
    end
    
    % 2. Bring in the data (d = light intensity, t = timestamps)
    data = load(dataFile);
    d = data.d; t = data.t; 
    
    % Hardcoded coefficients (E) and pathlength factor (DPF) 
    % Note: These should ideally be parsed from your .csv/txt files
    E = [0.0955, 0.4323;   
         0.2526, 0.1798];
    invE = inv(E); % Invert the matrix to solve for concentration later
    
    DPF = 6; 
    distance = SDS * DPF; % Total path light travels through the head
    
    % 3. Convert raw intensity to Optical Density (OD)
    % We use the average of the signal as a baseline
    baseline = mean(d, 1);
    OD = log10(repmat(baseline, size(d,1), 1) ./ d); 
    
    [nTime, nChan] = size(OD);
    half = nChan / 2; % NIRS usually has two wavelengths per location
    
    dHbO = zeros(nTime, half);
    dHbR = zeros(nTime, half);
    
    % 4. Apply the Beer-Lambert Law 
    % Loop through each source-detector pair to find HbO and HbR
    for i = 1:half
        % Pair up the two wavelengths for this channel
        temp_OD = [OD(:, i), OD(:, i+half)]'; 
        
        % The "magic" step: turning light absorption into hemoglobin concentration
        conc = invE * (temp_OD / distance);
        dHbO(:, i) = conc(1, :)';
        dHbR(:, i) = conc(2, :)';
    end
    
    % 5. Create the plots for the requested channels
    fig = figure;
    if ~isempty(plotChannelIdx)
        for j = 1:length(plotChannelIdx)
            ch = plotChannelIdx(j);
            subplot(length(plotChannelIdx), 1, j);
            
            % Red for Oxygenated, Blue for Deoxygenated
            plot(t, dHbO(:, ch), 'r', 'LineWidth', 1.5); hold on;
            plot(t, dHbR(:, ch), 'b', 'LineWidth', 1.5);
            
            title(['Concentration Changes - Channel ', num2str(ch)]);
            xlabel('Time (s)'); ylabel('\Delta Conc');
            legend('dHbO', 'dHbR');
            grid on; axis tight;
        end
    end
end
