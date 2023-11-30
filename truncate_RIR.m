function truncatedSignal = truncate_RIR(signal, threshold, start_offset)
    % Find the index of the first sample above the threshold
    startIndex = find(abs(signal) > threshold, 1, 'first');

    % If no sample is above the threshold, set startIndex to 1
    if isempty(startIndex)
        startIndex = 1;
    end

    % Truncate the signal from the startIndex
    truncatedSignal = signal(startIndex-start_offset:end);
end
