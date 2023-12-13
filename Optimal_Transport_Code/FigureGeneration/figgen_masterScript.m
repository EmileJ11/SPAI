%% SRIR Interpolation Examples: Figure Generation
% combines all the necessary stuff for important things
% Created by: Aaron Geldert
% Last modified: 5 Dec 2022

% The script produces figures:
% 1: room geometry w/src & rcv positions
% 2: room geometry and VS, 4 views
% 3: projection map of ground truth
% 4: projection map with OT between ground truth
% 5: projection map with NN between ground truth
% 6: rir plot of ground truth
% 7: rir plot of linear methods
% 8: rir plot of mapping methods

% %%%%%% 0) Initialize params   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

% metadata
numSetups = 12;                  % how many runs to do?

% GA Parameters
geomFile = 'cuboid.txt';   % which room geometry to use
maxIsOrder = 3;                 % image method maximum order
rCoeff = 0.707;                 % reflection coeff of all surfaces
interpDist = 2;                 % distance the rcv moves between P and Q
wallOffset = 0.5;               % gap btwn a wall and valid src/rcv position
global globalSurfaceCount

% OT Parameters
mapTol = 1.05;                  % cost (scaled to avoid equivalent dummy mappings)
xi = mapTol * interpDist.^2;    % dummy mapping cost
fs = 48e3;                      % sample rate Hz

% directory
parentFolder = '/Users/aarongeldert/Documents/MATLAB/THESIS/Figures/ThesisSVGs/';
parentFolder = [parentFolder erase(geomFile,'.txt') '_' num2str(interpDist) '/'];

%% RUN THIS SECTION REPEATEDLY!

for iterInd = 1:numSetups
    close all;
    disp(['%%%%%%%%%%%%%%%% SETUP ' num2str(iterInd) ' / ' num2str(numSetups) ':']);
    
% %%%%%% 1) Initialize IsmSimulator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
globalSurfaceCount = 1;
ISM = IsmSimulator(geomFile);
ISM.maxIsOrder = maxIsOrder;
ISM.rCoeff = rCoeff;
ISM.verbose = 0;

% Initialize source position
srcPos = NaN(size(1,3));
while ~ISM.isPointInsideGeometry(srcPos, wallOffset)
    srcPos = ISM.furthestPt .* rand(1,3); % find a random point inside
end
ISM.sourcePosition = srcPos;
ISM.simulateImageSources;

% Initialize receiver positions
% find 2 random points inside
rcvPos1 = NaN(size(1,3));
rcvPos2 = NaN(size(1,3));
while ~ISM.isPointInsideGeometry(rcvPos1, wallOffset) || ~ISM.isPointInsideGeometry(rcvPos2, wallOffset)
    rcvPos1 = ISM.furthestPt .* rand(1,3);
%     dir = randomPointOnSphere; % random in 3D
    dir = S2C([360*rand(1), 0, 1]); % random in 2D (horizontal)
    rcvPos2 = rcvPos1 + interpDist*dir;
end

% Initialize the POT solution

% Make params structure 
params.ISM = ISM;
params.xi = xi;
params.rcvPos1 = rcvPos1;
params.rcvPos2 = rcvPos2;
params.geomFile = geomFile;

rcvPosR = 0.5*rcvPos1 + 0.5*rcvPos2;

% %%%%%% Geometry plots 1-2     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
% ISM.receiverPosition = rcvPosR;
[hs, hr, hi] = ISM.plotGeometry(0); % no VS
hold on;
% add two receiver positions
hr1 = scatter3(rcvPos1(1), rcvPos1(2), rcvPos1(3), 80, 'g>', 'filled', 'MarkerEdgeColor', 'k');
line([rcvPos1(1) rcvPos1(1)],...
    [rcvPos1(2) rcvPos1(2)],...
    [0 rcvPos1(3)],'Color','black','LineStyle',':');
scatter3(rcvPos1(1), rcvPos1(2), 0,'k.');

hr2 = scatter3(rcvPos2(1), rcvPos2(2), rcvPos2(3), 80, 'g<', 'filled', 'MarkerEdgeColor', 'k');
line([rcvPos2(1) rcvPos2(1)],...
    [rcvPos2(2) rcvPos2(2)],...
    [0 rcvPos2(3)],'Color','black','LineStyle',':');
scatter3(rcvPos2(1), rcvPos2(2), 0,'k.');

hrR = scatter3(rcvPosR(1), rcvPosR(2), rcvPosR(3), 80, 'gd', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.25);
line([rcvPosR(1) rcvPosR(1)],...
    [rcvPosR(2) rcvPosR(2)],...
    [0 rcvPosR(3)],'Color','black','LineStyle',':');
scatter3(rcvPosR(1), rcvPosR(2), 0,'k.');

line([rcvPos1(1) rcvPos2(1)],...
    [rcvPos1(2) rcvPos2(2)],...
    [rcvPos1(3) rcvPos2(3)],'Color','#77AC30','LineStyle','--');

legend([hs hr1 hr2 hrR],'Source','Receiver at $\mathbf{s}_{\mathcal{P}}$',...
    'Receiver at $\mathbf{s}_{\mathcal{Q}}$','Receiver at $\mathbf{s}_{\mathcal{R}}$',...
    'Location','northwest');

% Virtual source + geometry plots
PCr = ISM.imageSourcesForRcvPos(rcvPosR);
figure(2); 
t = tiledlayout(2,2); nexttile;
[~, hr, hi] = ISM.plotGeometry(1);
view([60 22.5]);
legend([hi, hr], 'VS','Receiver','Location','northeast');

nexttile;
ISM.plotGeometry(1);
view([90 90]);

nexttile;
ISM.plotGeometry(1);
view([0 0]);

nexttile;
ISM.plotGeometry(1);
view([90 0]);

t.TileSpacing = 'compact';
t.Padding = 'compact';

% %%%%%% Projection plots 3-5   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figgenFunc_VSmaps(params, [3 4 5]);

% %%%%%% RIR plots 6-8          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figgenFunc_RIRplots(params, [6 7 8]);

% SAVING
folderName = [erase(geomFile,'.txt') '_ex' num2str(iterInd)];
mkdir(parentFolder,folderName);
figure(1);
set(gcf,'Position',[1    53   365   300]);
saveas(gcf, [parentFolder folderName '/geometry.png']);
saveas(gcf, [parentFolder folderName '/geometry.svg']);
figure(2);
set(gcf,'Position',[227   901   820   500]);
saveas(gcf, [parentFolder folderName '/virtualSources.png']);
saveas(gcf, [parentFolder folderName '/virtualSources.svg']);
figure(3);
set(gcf,'Position',[227        1355         820         500]);
saveas(gcf, [parentFolder folderName '/projTrue.png']);
saveas(gcf, [parentFolder folderName '/projTrue.svg']);
figure(4);
set(gcf,'Position',[-633        1355         820         500]);
saveas(gcf, [parentFolder folderName '/projOT.png']);
saveas(gcf, [parentFolder folderName '/projOT.svg']);
figure(5);
set(gcf,'Position',[-633   901   820   500]);
saveas(gcf, [parentFolder folderName '/projNN.png']);
saveas(gcf, [parentFolder folderName '/projNN.svg']);
figure(6);
set(gcf,'Position',[40   426   700   300]);
saveas(gcf, [parentFolder folderName '/rirTrue.png']);
saveas(gcf, [parentFolder folderName '/rirTrue.svg']);
figure(7);
set(gcf,'Position',[741   427   700   300]);
saveas(gcf, [parentFolder folderName '/rirLins.png']);
saveas(gcf, [parentFolder folderName '/rirLins.svg']);
figure(8);
set(gcf,'Position',[741    53   700   300]);
saveas(gcf, [parentFolder folderName '/rirMaps.png']);
saveas(gcf, [parentFolder folderName '/rirMaps.svg']);
disp('Saved .png files');
end