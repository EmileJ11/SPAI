function outputSignal = rt(inputSignal, thresholdDecay, windowSize)
    % Input:
    %   inputSignal: Input signal to be truncated.
    %   thresholdDecay: The desired decay in dB (e.g., -60 for RT60).
    %   windowSize: Number of samples to consider for comparison.

    % Output:
    %   outputSignal: Truncated output signal.

    % Find the peak value of the input signal
    peakValue = max(abs(inputSignal))


    % Initialize output signal
    outputSignal = zeros(size(inputSignal));

    % Iterate through the input signal
    for i = 20000:length(inputSignal)
        % Check if the current sample and the previous samples are below the threshold
        if i >= windowSize
            windowedSamples = inputSignal(i-windowSize+1 : i);
            %disp('max')
            %disp(max(abs(windowedSamples)))
            %disp('A/A db')
            %disp(20*log10(peakValue/max(abs(windowedSamples))))
            if 20*log10(peakValue/rms(windowedSamples)) >= thresholdDecay
                % If the maximum value in the window drops below the threshold, truncate the signal
                outputSignal = inputSignal(1:i);
                break;
            end
        end

        % Assign the current sample to the output signal
        outputSignal(i) = inputSignal(i);
    end
end
