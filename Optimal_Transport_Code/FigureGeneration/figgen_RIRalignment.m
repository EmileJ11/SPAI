%% Mono RIR interpolation comparisons
% Created by: Aaron Geldert
% Last modified: 28 Nov 2022

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
rCoeff = 0.707;

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
    rcvDir = S2C([randAz, 0, 1]);
    rcvPos2 = rcvPos1 + interpDist*rcvDir;
end

% Define the interpolation problem:
k = 0.5;

rcvPos2 = rcvPos1 + interpDist*rcvDir./vecnorm(rcvDir);
rcvPosk = (1-k)*rcvPos1 + k*rcvPos2;
PC1 = ISM.imageSourcesForRcvPos(rcvPos1); % Pos 1 
% figure(9); ISM.plotGeometry;
PC2 = ISM.imageSourcesForRcvPos(rcvPos2); % Pos 2
% figure(11); ISM.plotGeometry;
PCgt = ISM.imageSourcesForRcvPos(rcvPosk); % Pos GT
% figure(10); ISM.plotGeometry;

PCli = combinePCs(PC1, PC2, k); % linear 
PCal = combineAlignedPCs(PC1, PC2, k); % aligned linear

OT = PCOT(PC1, PC2, xi);
PCnn = OT.greedyInterp(k); % nearest neighbor
PCot = OT.interpPC(k); % partial OT

% set up mono RIRs
sf = hann(round(0.002 * fs));
[rir(:,1), tvec] = monoIRfromPC(PC1, 0.2, fs);
[rir(:,2), ~] = monoIRfromPC(PC2, 0.2, fs);
[rir(:,3), ~] = monoIRfromPC(PCgt, 0.2, fs);
[rir(:,4), ~] = monoIRfromPC(PCot, 0.2, fs);
[rir(:,5), ~] = monoIRfromPC(PCnn, 0.2, fs);
[rir(:,6), ~] = monoIRfromPC(PCal, 0.2, fs);
[rir(:,7), ~] = monoIRfromPC(PCli, 0.2, fs);

edc = zeros(size(rir));
for ii = 1:size(rir,2)
    edc(:,ii) = conv(rir(:,ii),sf,'same');
end

plotTime = 0.04; % ms of plot desired
dbGap = 35;
dbBln = 30; linBln = 10.^(dbBln/20);
dbOffset = 3;
dbTickGap = -20;
axisLim1 = [0 1e3*plotTime -2*dbGap-dbBln 0];
stemWidth = 1.5;

dbRir = db(rir) + dbOffset;
dbRir(dbRir < -dbBln) = -Inf;
dbEdc = db(edc+eps) + dbOffset; % the eps ensures that there is no -Inf
dbEdc(dbEdc < -dbBln) = -dbBln;

greyColor = [.75 .75 .75];
tvecMs = 1e3*tvec;

%% Figure 1: Ground truth P, R, Q

figure(1); clf; h1 = axes;
set(gcf, 'Position', [100 520 700 300]);
set(h1, 'Units', 'Normalized'); % First change to normalized units.
set(h1, 'OuterPosition', [.02, .1, .99, .85]);

% 1: Position P
stem(1e3*tvec, dbRir(:,1),'filled','k',...
    'BaseValue', -dbBln, 'MarkerSize', 3, 'LineWidth', stemWidth);
hold on;
hFill = fill([0, tvecMs, max(tvecMs)],...
    [-dbBln-eps; dbEdc(:,1); -dbBln-eps],...
    greyColor, 'FaceAlpha', 0.15, 'EdgeAlpha',0.7);
axis(axisLim1); yticklabels({}); xticklabels({});
set(h1,'YTick',[]);

% 2: Position R (ground truth)
h2 = axes('position',h1.Position);
set(h2,'Color', 'none');
hold on;
stem(tvecMs, dbRir(:,3)-dbGap,'filled','k',...
    'BaseValue', -dbBln-dbGap, 'MarkerSize', 3, 'LineWidth', stemWidth);
hFill = fill([0, tvecMs, max(tvecMs)],...
    [-dbBln-dbGap-eps; dbEdc(:,3)-dbGap; -dbBln-dbGap-eps],...
    greyColor, 'FaceAlpha', 0.15, 'EdgeAlpha',0.7);

axis(axisLim1); yticklabels({}); xticklabels({}); set(h2,'YTick',[]);

% 3: Position Q
h3 = axes('position',h1.Position);
set(h3,'Color', 'none');
hold on;
stem(tvecMs, dbRir(:,2)-2*dbGap,'filled','k',...
    'BaseValue', -dbBln-2*dbGap, 'MarkerSize', 3, 'LineWidth', stemWidth);
hFill = fill([0, tvecMs, max(tvecMs)],...
    [-dbBln-2*dbGap-eps; dbEdc(:,2)-2*dbGap; -dbBln-2*dbGap-eps],...
    greyColor, 'FaceAlpha', 0.15, 'EdgeAlpha',0.7);
hold on;
yticks([-2*dbGap+dbTickGap -2*dbGap -dbGap+dbTickGap -dbGap dbTickGap 0]); 
yticklabels({num2str(dbTickGap), '0', num2str(dbTickGap), '0', num2str(dbTickGap), '0'});
axis(axisLim1); 
% xticklabels({}); 
ylabel('$h(t)$ [dB]','FontSize',14);

title('')
xlabel('$t$ [ms]','FontSize',14);

% labels done manually here
textPosX = 0.18;
annotation('textbox', [textPosX, 0.85, 0, 0], 'string', '${\mathcal{P}}$',...
    'FontSize',12);
annotation('textbox', [textPosX, 0.62, 0, 0], 'string', '${\mathcal{R}}$',...
    'FontSize',12);
annotation('textbox', [textPosX, 0.38, 0, 0], 'string', '${\mathcal{Q}}$',...
    'FontSize',12);

saveas(gcf, '/Users/aarongeldert/Documents/MATLAB/THESIS/Figures/ThesisNewFigures/trapezoidal_rirAlign_orig.png');

%% Figure 2: Linear Combination methods
figure(2); clf; h1 = axes;
set(gcf, 'Position', [100 520 700 300]);
set(h1, 'Units', 'Normalized'); % First change to normalized units.
set(h1, 'OuterPosition', [.02, .1, .99, .85]);

% 1: Ground Truth
hFill = fill([0, tvecMs, max(tvecMs)],...
    [-dbBln-eps; dbEdc(:,3); -dbBln-eps],...
    greyColor, 'FaceAlpha', 0.15, 'EdgeAlpha',0.7);
hold on;
stem(1e3*tvec, dbRir(:,3),'filled','k',...
    'BaseValue', -dbBln, 'MarkerSize', 3, 'LineWidth', stemWidth);
axis(axisLim1); yticklabels({}); xticklabels({});
set(h1,'YTick',[]);

% 2: Aligned Linear Combination
h2 = axes('position',h1.Position);
set(h2,'Color', 'none');
hold on;
hFill = fill([0, tvecMs, max(tvecMs)],...
    [-dbBln-dbGap-eps; dbEdc(:,6)-dbGap; -dbBln-dbGap-eps],...
    [0.8500 0.3250 0.0980], 'FaceAlpha', 0.15, 'EdgeAlpha',0.7);
hold on;
stem(tvecMs, dbRir(:,6)-dbGap,'filled','Color','#D95319',...
    'BaseValue', -dbBln-dbGap, 'MarkerSize', 3, 'LineWidth', stemWidth);
axis(axisLim1); yticklabels({}); xticklabels({}); set(h2,'YTick',[]);

% 3: Linear Combination
h3 = axes('position',h1.Position);
set(h3,'Color', 'none');
hold on;
hFill = fill([0, tvecMs, max(tvecMs)],...
    [-dbBln-2*dbGap-eps; dbEdc(:,7)-2*dbGap; -dbBln-2*dbGap-eps],...
    [0.6350 0.0780 0.1840], 'FaceAlpha', 0.15, 'EdgeAlpha',0.7);
hold on;
stem(tvecMs, dbRir(:,7)-2*dbGap,'filled','Color','#A2142F',...
    'BaseValue', -dbBln-2*dbGap, 'MarkerSize', 3, 'LineWidth', stemWidth);
yticks([-2*dbGap+dbTickGap -2*dbGap -dbGap+dbTickGap -dbGap dbTickGap 0]); 
yticklabels({num2str(dbTickGap), '0', num2str(dbTickGap), '0', num2str(dbTickGap), '0'});
axis(axisLim1); 
% xticklabels({}); 
ylabel('$h(t)$ [dB]','FontSize',14);

title('')
xlabel('$t$ [ms]','FontSize',14);

% labels done manually here
textPosX = 0.18;
annotation('textbox', [textPosX, 0.85, 0, 0], 'string', '${\mathcal{R}}$',...
    'FontSize',12);
annotation('textbox', [textPosX, 0.62, 0, 0], 'string', '${\widehat{\mathcal{R}}}$, aligned linear',...
    'FontSize',12);
annotation('textbox', [textPosX, 0.38, 0, 0], 'string', '${\widehat{\mathcal{R}}}$, linear',...
    'FontSize',12);

saveas(gcf, '/Users/aarongeldert/Documents/MATLAB/THESIS/Figures/ThesisNewFigures/trapezoidal_rirAlign_lins.png');


%% Figure 3: Mapping methods
figure(3); clf; h1 = axes;
set(gcf, 'Position', [100 520 700 300]);
set(h1, 'Units', 'Normalized'); % First change to normalized units.
set(h1, 'OuterPosition', [.02, .1, .99, .85]);

% 1: Ground Truth
hFill = fill([0, tvecMs, max(tvecMs)],...
    [-dbBln-eps; dbEdc(:,3); -dbBln-eps],...
    greyColor, 'FaceAlpha', 0.15, 'EdgeAlpha',0.7);
hold on;
stem(1e3*tvec, dbRir(:,3),'filled','k',...
    'BaseValue', -dbBln, 'MarkerSize', 3, 'LineWidth', stemWidth);
axis(axisLim1); yticklabels({}); xticklabels({});
set(h1,'YTick',[]);

% 2: Partial Optimal Transport
h2 = axes('position',h1.Position);
set(h2,'Color', 'none');
hold on;
hFill = fill([0, tvecMs, max(tvecMs)],...
    [-dbBln-dbGap-eps; dbEdc(:,4)-dbGap; -dbBln-dbGap-eps],...
    [0 0.4470 0.7410], 'FaceAlpha', 0.15, 'EdgeAlpha',0.7);
hold on;
stem(tvecMs, dbRir(:,4)-dbGap,'filled','Color','#0072BD',...
    'BaseValue', -dbBln-dbGap, 'MarkerSize', 3, 'LineWidth', stemWidth);

axis(axisLim1); yticklabels({}); xticklabels({}); set(h2,'YTick',[]);

% 3: Greedy Assignment Mapping
h3 = axes('position',h1.Position);
set(h3,'Color', 'none');
hold on;
hFill = fill([0, tvecMs, max(tvecMs)],...
    [-dbBln-2*dbGap-eps; dbEdc(:,5)-2*dbGap; -dbBln-2*dbGap-eps],...
    [0.4660 0.6740 0.1880], 'FaceAlpha', 0.15, 'EdgeAlpha',0.7);
hold on;
stem(tvecMs, dbRir(:,5)-2*dbGap,'filled','Color','#77AC30',...
    'BaseValue', -dbBln-2*dbGap, 'MarkerSize', 3, 'LineWidth', stemWidth);
yticks([-2*dbGap+dbTickGap -2*dbGap -dbGap+dbTickGap -dbGap dbTickGap 0]); 
yticklabels({num2str(dbTickGap), '0', num2str(dbTickGap), '0', num2str(dbTickGap), '0'});
axis(axisLim1); 
% xticklabels({}); 
ylabel('$h(t)$ [dB]','FontSize',14);

title('')
xlabel('$t$ [ms]','FontSize',14);

% labels done manually here
textPosX = 0.18;
annotation('textbox', [textPosX, 0.85, 0, 0], 'string', '${\mathcal{R}}$',...
    'FontSize',12);
annotation('textbox', [textPosX, 0.62, 0, 0], 'string', '${\widehat{\mathcal{R}}}$, partial OT',...
    'FontSize',12);
annotation('textbox', [textPosX, 0.38, 0, 0], 'string', '${\widehat{\mathcal{R}}}$, greedy map',...
    'FontSize',12);

saveas(gcf, '/Users/aarongeldert/Documents/MATLAB/THESIS/Figures/ThesisNewFigures/trapezoidal_rirAlign_maps.png');


%%
% figure(2); clf; h3 = axes;
% set(gcf, 'Position', [100 100 600 310]);
% set(h3, 'Units', 'Normalized'); % First change to normalized units.
% set(h3, 'OuterPosition', [.02, .1, .99, .85]);
% 
% % 3: Position R
% stem(tvecMs, dbRir(:,3),'filled','k',...
%     'BaseValue', -dbBln, 'MarkerSize', 3, 'LineWidth', stemWidth);
% hold on;
% % plot(tvecMs, db(sir(:,1))+5,'k');
% axis(axisLim2); yticklabels({}); xticklabels({});
% set(h3,'YTick',[]);
% 
% rectangle('Position',[0 -dbBln 1e3*plotTime dbBln],'FaceColor',[0 0 0 .05],'EdgeColor','none');
% hold on;
% 
% % 4: Estimate POT
% h4 = axes('position',h3.Position);
% set(h4,'Color', 'none');
% hold on;
% stem(tvecMs, dbRir(:,3)-dbGap,'Color', greyColor,...
%     'BaseValue', -dbBln-dbGap, 'MarkerSize', 3, 'LineWidth', 1);
% hold on;
% stem(tvecMs, dbRir(:,4)-dbGap,'filled','Color','#0072BD',...
%     'BaseValue', -dbBln-dbGap, 'MarkerSize', 3, 'LineWidth', stemWidth);
% % plot(tvecMs, db(sir(:,2))-dbGap+5,'k');
% axis(axisLim2); xticklabels({}); yticklabels({}); set(h4,'YTick',[]);
% 
% % 5: Estimate NN
% h3 = axes('position',h4.Position);
% set(h3,'Color', 'none');
% hold on;
% stem(tvecMs, dbRir(:,3)-2*dbGap,'Color', greyColor,...
%     'BaseValue', -dbBln-2*dbGap, 'MarkerSize', 3, 'LineWidth', 1);
% hold on;
% stem(tvecMs, dbRir(:,5)-2*dbGap,'filled','Color','#77AC30',...
%     'BaseValue', -dbBln-2*dbGap, 'MarkerSize', 3, 'LineWidth', stemWidth);
% yticks([-3*dbGap+dbTickGap -3*dbGap -2*dbGap+dbTickGap -2*dbGap -dbGap+dbTickGap -dbGap dbTickGap 0]); 
% yticklabels({num2str(dbTickGap), '0', num2str(dbTickGap), '0',...
%     num2str(dbTickGap), '0',num2str(dbTickGap), '0'});
% axis(axisLim2); xticklabels({}); 
% ylabel('$h$ (dB)','FontSize',14);
% 
% % 6: Estimate Linear
% h6 = axes('position',h4.Position);
% set(h6,'Color', 'none');
% hold on;
% stem(tvecMs, dbRir(:,3)-3*dbGap,'Color', greyColor,...
%     'BaseValue', -dbBln-3*dbGap, 'MarkerSize', 3, 'LineWidth', 1);
% hold on;
% stem(tvecMs, dbRir(:,6)-3*dbGap,'filled','Color','#A2142F',...
%     'BaseValue', -dbBln-3*dbGap, 'MarkerSize', 3, 'LineWidth', stemWidth);
% hold on;
% % plot(tvecMs, db(sir(:,3))-2*dbGap+5,'k');
% axis(axisLim2); yticklabels({});
% set(h6,'YTick',[]);
% 
% annotation('textbox', [textPosX, 0.86, 0, 0], 'string', 'c)',...
%     'FontSize',12);
% annotation('textbox', [textPosX, 0.71, 0, 0], 'string', 'd)',...
%     'FontSize',12);
% annotation('textbox', [textPosX, 0.54, 0, 0], 'string', 'e)',...
%     'FontSize',12);
% annotation('textbox', [textPosX, 0.36, 0, 0], 'string', 'f)',...
%     'FontSize',12);
% xlabel('$t$ (ms)','FontSize',14);
