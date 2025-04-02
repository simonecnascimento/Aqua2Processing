function plotFFT(signal, Fs)
    % Function to plot the FFT of a signal
    % Inputs:
    %   signal - The time-domain signal (vector)
    %   Fs - Sampling frequency (Hz)

    %1.Identify Signal and Noise Components – The FFT helps you analyze the frequency content of your signal, making it easier to determine if there are dominant frequencies and how much noise is present.
    %2.Optimize Filtering Parameters – Based on the frequency spectrum, you can choose an appropriate Savitzky-Golay filter (sgolayfilt) frame size and polynomial order to remove noise effectively.
    %3.Ensure Proper SNR Calculation – By understanding the frequency content beforehand, you can better define what constitutes "signal" and "noise" in your analysis.
    
    % Compute the FFT
    L = length(signal);  % Length of signal
    Y = fft(signal);     % Compute FFT
    P2 = abs(Y/L);       % Normalize magnitude
    P1 = P2(1:L/2+1);    % Take only first half (positive frequencies)
    P1(2:end-1) = 2*P1(2:end-1);  % Adjust magnitude for single-sided spectrum
    
    % Frequency axis
    f = Fs * (0:(L/2)) / L;
    
    % Plot the frequency spectrum
    figure;
    plot(f, P1, 'b', 'LineWidth', 1.5);
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    title('Frequency Spectrum');
    grid on;
end
