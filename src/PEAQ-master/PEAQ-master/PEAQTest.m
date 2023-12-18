%%
clear; clc;

ref = 'Measured_speaker2.wav';
test = 'Auralized_speaker2.wav';
PQevalAudio(ref, test);

% Load the saved MOVC structure
load MOVC.mat;



%% Plot Modulation Difference (ModDiff) for each frame
FrameNumber = 1:length(MOVC.MDiff.Mt1B);
figure;
plot(FrameNumber, MOVC.MDiff.Mt1B, 'b', 'LineWidth', 2, 'DisplayName', 'Mt1B');
hold on;
plot(FrameNumber, MOVC.MDiff.Mt2B, 'r', 'LineWidth', 2, 'DisplayName', 'Mt2B');
xlabel('Frame Number');
ylabel('Modulation Difference');
legend('Mt1B', 'Mt2B');
title('Modulation Difference Over Frames');
grid on;
hold off;

%% Plot Non-Stationary-to-Stationary Ratio (NMR) for each frame
FrameNumber = 1:length(MOVC.MDiff.Mt1B);
figure;
plot(FrameNumber, MOVC.NMR.NMRavg, 'g', 'LineWidth', 2, 'DisplayName', 'NMRavg');
hold on;
plot(FrameNumber, MOVC.NMR.NMRmax, 'm', 'LineWidth', 2, 'DisplayName', 'NMRmax');
xlabel('Frame Number');
ylabel('NMR');
legend('NMRavg', 'NMRmax');
title('Non-Stationary-to-Stationary Ratio (NMR) Over Frames');
grid on;
hold off;

%% Plot Bandwidth (BW) for Reference and Test Signals for each frame
FrameNumber = 1:length(MOVC.MDiff.Mt1B);
figure;
plot(FrameNumber, MOVC.BW.BWRef, 'c', 'LineWidth', 2, 'DisplayName', 'BWRef');
hold on;
plot(FrameNumber, MOVC.BW.BWTest, 'y', 'LineWidth', 2, 'DisplayName', 'BWTest');
xlabel('Frame Number');
ylabel('Bandwidth');
legend('BWRef', 'BWTest');
title('Bandwidth (BW) Over Frames');
grid on;
hold off;

%% Plot Loudness for Reference and Test Signals for each frame
FrameNumber = 1:length(MOVC.MDiff.Mt1B);
figure;
plot(FrameNumber, MOVC.Loud.NRef, 'b', 'LineWidth', 2, 'DisplayName', 'NRef');
hold on;
plot(FrameNumber, MOVC.Loud.NTest, 'r', 'LineWidth', 2, 'DisplayName', 'NTest');
xlabel('Frame Number');
ylabel('Loudness');
legend('NRef', 'NTest');
title('Loudness Over Frames');
grid on;
hold off;

%% Plot Perceptual Distortion (PD) for each frame
FrameNumber = 1:length(MOVC.MDiff.Mt1B);
figure;
plot(FrameNumber, MOVC.PD.Pc, 'g', 'LineWidth', 2, 'DisplayName', 'Pc');
hold on;
plot(FrameNumber, MOVC.PD.Qc, 'm', 'LineWidth', 2, 'DisplayName', 'Qc');
xlabel('Frame Number');
ylabel('Perceptual Distortion');
legend('Pc', 'Qc');
title('Perceptual Distortion (PD) Over Frames');
grid on;
hold off;

%% Plot Equivalent Harmonic Signals (EHS) for each frame
FrameNumber = 1:length(MOVC.MDiff.Mt1B);
figure;
plot(FrameNumber, MOVC.EHS.EHS, 'k', 'LineWidth', 2, 'DisplayName', 'EHS');
xlabel('Frame Number');
ylabel('EHS');
legend('EHS');
title('Equivalent Harmonic Signals (EHS) Over Frames');
grid on;
hold off;

%% Plotting single-value metrics 

metrics = [-11.8855, -1.331]; % Example values
metricNames = {'Total NMRB', 'Objective Difference Grade'};

figure;
bar(metrics);
set(gca, 'xticklabel', metricNames);
ylabel('Value');
title('Single-Value Metrics');

%% Windowed modulation difference

FrameNumber = 1:length(MOVC.MDiff.Mt1B);
figure;
plot(FrameNumber, MOVC.MDiff.Wt, 'DisplayName', 'WinModDiff1B'); % Replace with actual frame-based data
hold on;
xlabel('Frame Number');
ylabel('Value');
legend('show');
title('Windowed modulation difference');
grid on;
hold off;
