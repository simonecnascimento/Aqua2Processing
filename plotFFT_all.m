function plotFFT_all(dFF_all, Fs)
    % Function to compute and plot the histogram of dominant frequencies for all cells
    % Inputs:
    %   dFF_all - Matrix of time-domain signals for all cells (each row is a signal)
    %   Fs - Sampling frequency (Hz)
    
    % Loop over all cells to compute FFT magnitudes
    numCells = size(dFF_all, 1);
    dominantFrequencies = zeros(numCells, 1);

    for w = 1:numCells
        signal = dFF_all(w, :);
        
        % Compute the FFT
        L = length(signal);  % Length of signal
        Y = fft(signal);     % Compute FFT
        P2 = abs(Y/L);       % Normalize magnitude
        P1 = P2(1:L/2+1);    % Take only first half (positive frequencies)
        P1(2:end-1) = 2*P1(2:end-1);  % Adjust magnitude for single-sided spectrum
        
        % Frequency axis
        f = Fs * (0:(L/2)) / L;

        % Find the frequency corresponding to the peak in the FFT magnitude
        [~, maxIdx] = max(P1);  % Index of maximum magnitude
        dominantFrequencies(w) = P1(maxIdx);  % Dominant frequency for this cell

    end

    % Plot histogram of the dominant frequencies
    figure;
    histogram(dominantFrequencies, 'BinMethod', 'auto', 'FaceColor', 'b', 'EdgeColor', 'k');
    xlabel('Dominant Frequency (Hz)');
    ylabel('Number of Cells');
    title('Histogram of Dominant Frequencies Across All Cells');
    grid on;
 
end