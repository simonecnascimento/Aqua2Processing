function plotFilteredSignal(signal, best_filtered_signal, best_frame_size, numTimePoints, i)
    % Set figure window to normal
    set(0, 'DefaultFigureWindowStyle', 'normal');
    
    figure;
    subplot(3,1,1);
    plot(signal, 'k', 'LineWidth', 1); hold on;
    plot(best_filtered_signal, 'r', 'LineWidth', 1.5);
    title(sprintf('Cell %d - Best Frame Size: %d', i, best_frame_size));
    legend('Original Signal', 'Filtered Signal');
    xlabel('Time'); ylabel('dF/F');

    subplot(3,1,2);
    plot(signal - best_filtered_signal, 'b');
    title('Noise (Original - Filtered)');
    xlabel('Time'); ylabel('Residual');

    subplot(3,1,3);
    fft_signal = abs(fft(signal)); % Compute FFT
    freq = linspace(0, 1, numTimePoints/2); % Frequency axis
    plot(freq, fft_signal(1:numTimePoints/2), 'g');
    title('FFT of Original Signal');
    xlabel('Frequency'); ylabel('Magnitude');

    pause(0.5); % Pause to view each plot
end
