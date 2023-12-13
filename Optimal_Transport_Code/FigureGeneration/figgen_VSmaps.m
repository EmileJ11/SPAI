%% Projections of virtual source interpolations
% Created by: Aaron Geldert
% Last modified: 28 Nov 2022
% better to use figgenFunc_VSmaps.m now!

clear; close all; clc; 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

% MODEL PARAMETERS
geomFile = 'trapezoidal.txt';
maxIsOrder = 3;
mapTol = 1.05;
interpDist = 2;
xi = mapTol * interpDist.^2;
wallOffset = 0.2;
rCoeff = 0.7071;

% setup the Simulator
global globalSurfaceCount
globalSurfaceCount = 1;
ISM = IsmSimulator(geomFile);
ISM.maxIsOrder = maxIsOrder;
ISM.rCoeff = rCoeff;
fs = ISM.fs;

% Source position
srcPos = NaN(size(1,3));
while ~ISM.isPointInsideGeometry(srcPos, wallOffset)
    srcPos = ISM.furthestPt .* rand(1,3); % find a random point inside
end

% now we have a valid path of interpolation inside the room geometry
ISM.sourcePosition = srcPos;
ISM.simulateImageSources;

% Receiver positions
% find 2 random points inside
rcvPos1 = NaN(size(1,3));
rcvPos2 = NaN(size(1,3));
while ~ISM.isPointInsideGeometry(rcvPos1, wallOffset) || ~ISM.isPointInsideGeometry(rcvPos2, wallOffset)
    rcvPos1 = ISM.furthestPt .* rand(1,3);
%     dir = randomPointOnSphere; % random in 3D
    randAz = 360*rand(1);
    dir = S2C([randAz, 0, 1]);
    rcvPos2 = rcvPos1 + interpDist*dir;
end

% Setup the Map Plot
figure(1);
aziGrid = [0 30  60 90 120 150 180 -30 -60  -120 -150 -90 179.9]'/180*pi;
eleGrid = [-60 -30 0 30 60 ]'/180*pi;
linegray = 0.75;
% subplot(6, 1, 1:3);
for k=1:length(aziGrid)
    ele=linspace(-pi/2, pi/2, 50);
    azi=ones(size(ele))*aziGrid(k);
    xyGrid = hammerAidhofProjection([azi', ele']);
    plot(xyGrid(:, 1), xyGrid(:, 2),'--', 'color', [1 1 1]*linegray);
    hold on
    text(xyGrid(25, 1),0, num2str(round(180/pi*aziGrid(k))), ...
        'fontsize', 9, 'color', [1 1 1]*linegray );
end
for k=1:length(eleGrid)
    azi=linspace(-pi, pi-.01, 50);
    ele=ones(size(azi))*eleGrid(k);
    xyGrid = hammerAidhofProjection([azi', ele']);
    plot(xyGrid(:, 1), xyGrid(:, 2),'--', 'color', [1 1 1]*linegray);
    hold on
    text(0, xyGrid(25, 2), num2str(180/pi*eleGrid(k)), ...
        'fontsize', 9, 'color', [1 1 1]*linegray );
end

% define smooth interpolation param values for graphics
numk = 251;
kVec = linspace(0, 1, numk);
lineScale = 4;

for kInd = 2:(numk-1)
    k = kVec(kInd);
    rcvPos = (1-k)*rcvPos1 + k*rcvPos2;
    PCk = ISM.imageSourcesForRcvPos(rcvPos);
    PCksph = deg2rad(C2S(PCk.pos));
    PCkhap = hammerAidhofProjection(PCksph(:,1:2));
    hk = scatter(PCkhap(:,1), PCkhap(:,2), lineScale, [.9 .6 .95], 'filled');
end

%% define coarse interpolation param values for comparisons
kVec = [0 0.5 1];
numk = numel(kVec);
rcvPos = NaN(numk,3);

% initialize PC structs
PCgt = struct('n',cell(1,numk),...
                'pos',cell(1,numk),...
                'mass',cell(1,numk));
for kInd = 1:numk
    k = kVec(kInd);
    rcvPos(kInd,:) = (1-k)*rcvPos1 + k*rcvPos2;
    PCgt(kInd) = ISM.imageSourcesForRcvPos(rcvPos(kInd,:));
end

% Generate 4 interpolated estimates
k = 0.5;
PCli = combinePCs(PCgt(1), PCgt(end), k); % linear 
PCal = combineAlignedPCs(PCgt(1), PCgt(end), k); % aligned linear

% OT = PCOT(PCgt(1), PCgt(end), xi);
OT = VSOT(PCgt(1), PCgt(end), interpDist, 0.05);
PCnn = OT.greedyInterp(k); % nearest neighbor
PCot = OT.interpPC(k); % partial OT

% make struct of PCs
PC1 = PCgt(1);
PC2 = PCgt(end);
% PCk = PCot;

% Plot VSs
scale = 25; 

PC1sph = deg2rad(C2S(PC1.pos));
PC1hap = hammerAidhofProjection(PC1sph(:,1:2));
h1 = scatter(PC1hap(:,1), PC1hap(:,2), lineScale+1.5+db(1+scale*PC1.mass).^1.5, [0 0 1], 'filled', 'MarkerFaceAlpha', .6);

PC2sph = deg2rad(C2S(PC2.pos));
PC2hap = hammerAidhofProjection(PC2sph(:,1:2));
h2 = scatter(PC2hap(:,1), PC2hap(:,2), lineScale+1.5+db(1+scale*PC2.mass).^1.5, [1 0 0], 'filled', 'MarkerFaceAlpha',.6);
% 
% PCksph = deg2rad(C2S(PCk.pos));
% PCkhap = hammerAidhofProjection(PCksph(:,1:2));
% scatter(PCkhap(:,1), PCkhap(:,2), 1+db(lineScale+scale*PCk.mass).^1.5, [.5 0 .5], 'filled');

% title('Projection map of Virtual Sources, receiver perspective');
xlabel('Azimuth ($^{\circ}$)'); xticks([]);
ylabel('Elevation ($^{\circ}$)'); yticks([]);
set(gcf, 'Position', [100 100 900 450]);

legend([h1 hk h2],'$\mathcal{P}$','$\int_{0}^{1} \mathcal{R(\kappa)} \mathrm{d}\kappa$','$\mathcal{Q}$', 'Location','southwest');
title('\textbf{Trapezoidal room}');
