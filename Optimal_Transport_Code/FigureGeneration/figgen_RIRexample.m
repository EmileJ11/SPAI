%% RIR Plot and Echogram

clear; clc; close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
[x, fs] = audioread('reference_arni_rir.wav');
x = x(:,1);
t = 0:1:(numel(x)-1);
t = t./fs;

% RIR
subplot(121);
plot(t,x,'k-');
axis([0 0.25 -1 1]);
xlabel('Time (sec)');
ylabel('Amplitude');

% Echogram
en = 10*log10((x+eps).^2);
[B,A] = butter(8, 0.1);
en = filtfilt(B,A,en);

subplot(122);
plot(t,en,'k-');
axis([0 0.25 -40 0]);
xlabel('Time (sec)');
ylabel('dB')

set(gcf,'Position',[100 100 800 220]);
saveas(gcf, '/Users/aarongeldert/Documents/MATLAB/THESIS/Figures/ThesisNewFigures/rir_mono.png');
