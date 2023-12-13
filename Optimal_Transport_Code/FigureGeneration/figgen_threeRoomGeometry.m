%% Generate 3 side-by-side images of the test rooms
% Created by: Aaron Geldert
% Last modified: 1 Dec 2022

clear; close all; clc;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

global globalSurfaceCount
globalSurfaceCount = 1;
axisLim = [-1 9 -1 11 0 5.5];
camView = [-30 22.5];

t = tiledlayout(1,3);
nexttile;
ISM = IsmSimulator('cuboid.txt');
ISM.plotGeometry;
title('\textbf{Cuboid}');
axis(axisLim);
view(camView);
grid on;

nexttile;
ISM = IsmSimulator('canted.txt');
ISM.plotGeometry;
title('\textbf{Canted}');
axis(axisLim);
view(camView);
grid on;
xlabel(''); ylabel(''); zlabel('');

nexttile;
ISM = IsmSimulator('trapezoidal.txt');
ISM.plotGeometry;
title('\textbf{Trapezoidal}');
axis(axisLim);
view(camView);
grid on;
xlabel(''); ylabel(''); zlabel('');

t.TileSpacing = 'compact';
t.Padding = 'compact';
set(gcf,'Position',[100 100 750 250]);

