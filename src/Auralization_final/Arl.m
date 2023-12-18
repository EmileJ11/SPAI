% Load RIRs (adjust the file path and variable name as needed)
load('rir_rt40.mat');
rirs = speakers_all_discrete_positions;

% Load HRTFs (adjust the file path and variable name as needed)
load('KU100_CDFC.mat');
hrtfs = hpcf.minPhase;


%% Convolve RIR with HRTF

numConvolvedSignals = 0;
maxSignalLength = size(hrtfs, 1) + size(rirs{1}, 1) - 1;
SumSignal = zeros(maxSignalLength, 1);


% Iterate through each RIR
for i = 1:numel(rirs)
    rir = rirs{i};

    % Assuming rir is a matrix
    for k = 1:size(rir, 2)
         rircolumn = rir(:, k);

        % Convolve with HRTF
        convolvedSignal = conv(hrtfs, rircolumn);

        % Normalize the convolved signal
        convolvedSignal = convolvedSignal / max(abs(convolvedSignal));

        % Zero-pad or trim the convolved signal to match SumSignal's size
        if length(convolvedSignal) < maxSignalLength
            % Zero-pad if convolvedSignal is shorter
            convolvedSignal(end:maxSignalLength) = 0;
        else
            % Trim if convolvedSignal is longer
            convolvedSignal = convolvedSignal(1:maxSignalLength);
        end

        % Update SumSignal and numConvolvedSignals
        numConvolvedSignals = numConvolvedSignals + 1;
        SumSignal = SumSignal + convolvedSignal;
    end
end

SumSignal = SumSignal/numConvolvedSignals;

%% Sound file

outputFileName = 'Binaural.wav';

% Specify the sample rate (assuming it is the same as the original signal)
sampleRate = 48000;  % Change this if your sample rate is different

% Save the omnidirectionalSignal as a .wav file
audiowrite(outputFileName, SumSignal, sampleRate);

%%






