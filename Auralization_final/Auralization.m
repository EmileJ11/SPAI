%% Load all of the data from the SRIR and the HRTF files first.
clear;
clc;

% Load HRTF data
load('hrir_final.mat'); 

% Load SRIR data
load('DOA_rt40.mat');    
load('pressure_rt40.mat');


% Initialize the BRIRs with zeros for faster calculation.
brir_left = cell(12, 1);
brir_right = cell(12, 1);

for i = 1:12
    brir_left{i} = zeros(size(P{1,1}));
    brir_right{i} = zeros(size(P{1,1}));
end

%% Making the BRIRs for left and right ear. 

speaker = 2; % Change speaker here!

for i = 1:size(DOA, 1) % Loop through rows
        currentElement = DOA{i, speaker}; % Change the number based on speaker
        [az, elv] = cart2sph(currentElement(:,1),currentElement(:,2),currentElement(:,3));
        azimuth = reshape(OnR, [], 1);
        elevation = reshape(OnL, [], 1);
        HRTF_positions = [azimuth(:), elevation(:)];
        query_points = [az(:), elv(:)];
        query_points_in_degrees(:, 1) = rad2deg(query_points(:, 1));
        query_points_in_degrees(:, 2) = rad2deg(query_points(:, 2));
        nearest_indices = knnsearch(HRTF_positions, query_points_in_degrees);
        nearest_hrtfs_left = hrir_l(nearest_indices);
        nearest_hrtfs_right = hrir_r(nearest_indices);
        srirs = P{i,speaker}; %Change the number based on speaker 

        temp_brir_left = zeros(size(srirs));
        temp_brir_right = zeros(size(srirs));

        for k = 1:size(srirs, 1)
            temp_brir_left(k, :) = conv(srirs(k, :), nearest_hrtfs_left(k, :));
            temp_brir_right(k, :) = conv(srirs(k, :), nearest_hrtfs_right(k, :));
        end

        brir_left{i} = temp_brir_left;
        brir_right{i} = temp_brir_right;
end

%% Final auralization

[audio_data, Fs] = audioread('OG_sound.wav');
[audio_data2, Fs2] = audioread('continuous_speaker1_ADAT8.wav');

% Determine the length of the shorter file (in samples)
length_audio2 = size(audio_data2, 1);

% Truncate the first audio file to the same length
audio_data_truncated = audio_data(1:length_audio2, :);

% Initialize convolved signals for faster convolution
total_convolved_left = zeros(size(audio_data_truncated, 1), 1);
total_convolved_right = zeros(size(audio_data_truncated, 1), 1);

for i = 1:length(brir_left)
    current_brir_left = brir_left{i};
    current_brir_right = brir_right{i};

    if size(audio_data_truncated, 2) == 2
        % Stereo
        convolved_left = conv(audio_data_truncated(:, 1), current_brir_left, 'same');
        convolved_right = conv(audio_data_truncated(:, 2), current_brir_right, 'same');
    else
        % Mono
        convolved_left = conv(audio_data_truncated, current_brir_left, 'same');
        convolved_right = conv(audio_data_truncated, current_brir_right, 'same');
    end

    % Sum up the convolved signals
    total_convolved_left = total_convolved_left + convolved_left;
    total_convolved_right = total_convolved_right + convolved_right;
end

% Combine the channels
convolved_audio = [total_convolved_left, total_convolved_right];

convolved_audio = convolved_audio / max(abs(convolved_audio(:)));

audiowrite('Auralized_audioforstep6_speaker2.wav', convolved_audio, Fs);

%% Load the correct azimuth and elevation values from the hrir_final.mat file.

% Assuming OnR and OnL are the azimuth and elevation data
% azimuth = reshape(OnR, [], 1);
% elevation = reshape(OnL, [], 1);

%% Create arrays for the both the loaded azimuths and elevations and use them for knnsearch function.

% HRTF_positions = [azimuth(:), elevation(:)];
% 
% query_points = [az(:), elv(:)];
% query_points_in_degrees(:, 1) = rad2deg(query_points(:, 1));
% query_points_in_degrees(:, 2) = rad2deg(query_points(:, 2));
% 
% % Find the nearest HRTFs for each query point
% nearest_indices = knnsearch(HRTF_positions, query_points_in_degrees);
% 
% nearest_hrtfs_left = hrir_l(nearest_indices);
% nearest_hrtfs_right = hrir_r(nearest_indices);

%% Convolve the pressure signals with the corresponding nearest hrtfs left and right

% srirs = P{12,2};
% 
% for i = 1:size(srirs, 1)
%     hrtf_index = nearest_indices(i);
%     brir_left(i, :) = conv(srirs(i, :), nearest_hrtfs_left(i, :));
%     brir_right(i, :) = conv(srirs(i, :), nearest_hrtfs_right(i, :));
% end



