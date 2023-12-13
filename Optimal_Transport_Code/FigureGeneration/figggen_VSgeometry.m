%% Thesis figures of ISM Simulator

clear; close all; clc; 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
% Prepare the Simulator
global globalSurfaceCount
globalSurfaceCount = 1;


%% Figure 1: Rectangular room geometry and VS (3D views)
clear; close all; clc;
global globalSurfaceCount
globalSurfaceCount = 1;

rCoeff = 0.7071;
ISM = IsmSimulator('cuboid.txt');
ISM.maxIsOrder = 3;
ISM.rCoeff = rCoeff;
fs = ISM.fs;
ISM.sourcePosition = [5 1.1 2.2];
ISM.simulateImageSources;
ISM.imageSourcesForRcvPos([3.2 2.6 1.8]);

vol.rectVol = prod(ISM.furthestPt);
vol.rt60 = 0.161 * vol.rectVol / (ISM.totalSurfaceArea * (1-rCoeff));
disp(vol);

figure(1); 
subplot(221);
ISM.plotGeometry(1);
view([0 0]);

subplot(222);
ISM.plotGeometry(1);
view([90 0]);

subplot(223);
ISM.plotGeometry(1);
view([0 90]);

subplot(224);
[hs, hr, hi] = ISM.plotGeometry(1);
view([30 30]);
legend([hs, hi, hr], 'Source','VS','Receiver','Location','northeast');


set(gcf,'Position',[100 100 900 500]);
saveas(gcf, '/Users/aarongeldert/Documents/MATLAB/THESIS/Figures/ThesisNewFigures/rectangularVS3.png');


%% Figure 2: Trapezoidal room geometry and VS (3D views)
clear; close all; clc;
global globalSurfaceCount
globalSurfaceCount = 1;

rCoeff = 0.7071;
ISM = IsmSimulator('trapezoidal.txt');
ISM.maxIsOrder = 3;
ISM.rCoeff = rCoeff;
fs = ISM.fs;

ISM.sourcePosition = [5 5.4 2.6];
ISM.simulateImageSources;
ISM.imageSourcesForRcvPos([3.2 2.6 1.8]);

figure(2); 
subplot(221);
ISM.plotGeometry;
view([0 0]);

subplot(222);
ISM.plotGeometry;
view([90 0]);

subplot(223);
ISM.plotGeometry;
view([90 90]);

subplot(224);
[hs, hr, hi] = ISM.plotGeometry;
legend([hs, hi, hr], 'Source','VS','Receiver','Location','northwest');
view([60 30]);

set(gcf,'Position',[100 100 900 500]);
% saveas(gcf, '/Users/aarongeldert/Documents/MATLAB/THESIS/Figures/ThesisNewFigures/trapezoidalVS3.png');







