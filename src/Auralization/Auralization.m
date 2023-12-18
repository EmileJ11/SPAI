%% Load all of the data from the SRIR and the HRTF files first.
clear;
clc;

% Load HRTF data
load('hrir_final.mat'); 

% Load SRIR data
load('DOA_rt40.mat');    
load('pressure_rt40.mat');

brir_left = zeros(size(P{1,1}));
brir_right = zeros(size(P{1,1}));

%% The full process of making binaural rirs check the commented sections below for individual code

for i = 1:size(DOA, 1) % Loop through rows
        currentElement = DOA{i, 1}; % Change the number based on speaker
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
        srirs = P{i,1}; %Change the number based on speaker 
        for k = 1:size(srirs, 1)
            brir_left(k, :) = conv(srirs(k, :), nearest_hrtfs_left(k, :));
            brir_right(k, :) = conv(srirs(k, :), nearest_hrtfs_right(k, :));
        end
end

%% Final auralization

[audio_data, Fs] = audioread('OG_sound.wav');
[audio_data2, Fs2] = audioread('continuous_speaker1_ADAT8.wav');

% Determine the length of the shorter file (in samples)
length_audio2 = size(audio_data2, 1);

% Truncate the first audio file to the same length
audio_data_truncated = audio_data(1:length_audio2, :);

if size(audio_data_truncated, 2) == 2
    convolved_left = conv(audio_data_truncated(:, 1), brir_left, 'same');
    convolved_right = conv(audio_data_truncated(:, 2), brir_right, 'same');
    convolved_audio = [convolved_left, convolved_right];
else
    % If mono, just duplicate the convolved signal for both channels
    convolved_signal = conv(audio_data_truncated, binaural_rir_left, 'same'); % Assuming mono signal should be convolved with left RIR
    convolved_audio = [convolved_signal, convolved_signal];
end

convolved_audio = convolved_audio / max(abs(convolved_audio(:)));

audiowrite('Auralized_audio_speaker2.wav', convolved_audio, Fs);




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



