% Make Plots for Comparisons
% Created by: Aaron Geldert
% Last modified: 24 Nov 2022

addpath(genpath('legendflex'));

%% 2 meter trapezoidal room comparison
% Load Data
load('/Users/aarongeldert/Documents/MATLAB/THESIS/New_ISM_Data/trapezoidal3_dist2_mapTol1.05/Errors/ERR_trapezoidal_order3_dist2.mat');

% Check Data (?)

axisLim = [0 1 0 0.4];
legendBuff = [-5 -12];

% Setup figure
% set(groot,'defaultAxesFontName','Verdana');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
figure(1); clf; 
subplot(131);

% lin vs pot in DB
[hMlic, hFlic] = iqrPlot(kVec, ERR.lic, 'lin'); hold on;
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');
Legend = {'Linear','Partial OT'};
[~,object_h,~,~] = legendflex([hFlic, hFotc],Legend,...
                            'anchor',{'ne','ne'},...  
                            'xscale',0.55,...
                            'buffer',legendBuff,... 
                            'fontSize',11);
axis(axisLim); yticks(0:0.1:0.5);
ylabel('$\mathcal{E}$','Interpreter','latex','FontSize',15);
xticks([0 0.25 0.5 0.75 1]);
hatchfill2(object_h(end-1), 'single','HatchAngle',60,'HatchColor',[0.6350 0.0780 0.1840],'HatchDensity',13); 
hatchfill2(object_h(end), 'single','HatchAngle',120,'HatchColor',[0 0.4470 0.7410],'HatchDensity',20);
object_h(end-1).FaceAlpha = 0.15;
object_h(end).FaceAlpha = 0.15;
hMlic.LineWidth = 1.5;
hMotc.LineWidth = 2;
xlabel('$\kappa$','Interpreter','latex','FontSize',16);
set(gca, 'YGrid', 'on', 'XGrid', 'off')
                        
% aligned vs pot
subplot(132);
[hMalc, hFalc] = iqrPlot(kVec, ERR.alc, 'al'); hold on;
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');
Legend = {'Aligned Linear','Partial OT'};
[~,object_h,~,~] = legendflex([hFalc, hFotc],Legend,...
                            'anchor',{'ne','ne'},...  
                            'xscale',0.55,...
                            'buffer',legendBuff,... 
                            'fontSize',11);
axis(axisLim); yticks(0:0.1:0.5);
xticks([0 0.25 0.5 0.75 1]);
object_h(end-1).FaceAlpha = 0.15;
object_h(end).FaceAlpha = 0.15;
hatchfill2(object_h(end-1), 'single','HatchAngle',60,'HatchColor',[0.8500 0.3250 0.098],'HatchDensity',13); 
hatchfill2(object_h(end), 'single','HatchAngle',120,'HatchColor',[0 0.4470 0.7410],'HatchDensity',20);
hMalc.LineWidth = 1.5;
hMotc.LineWidth = 2;
xlabel('$\kappa$','Interpreter','latex','FontSize',16);
set(gca, 'YGrid', 'on', 'XGrid', 'off')

% nn vs pot
subplot(133);
[hMnnc, hFnnc] = iqrPlot(kVec, ERR.nnc, 'nn'); hold on;
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');
Legend = {'Nearest Neighbor','Partial OT'};
[~,object_h,~,~] = legendflex([hFnnc, hFotc],Legend,...
                            'anchor',{'ne','ne'},...  
                            'xscale',0.55,...
                            'buffer',legendBuff,... 
                            'fontSize',11);
axis(axisLim); yticks(0:0.1:0.5);
xticks([0 0.25 0.5 0.75 1]);
hatchfill2(object_h(end-1), 'single','HatchAngle',60,'HatchColor',[0.4660 0.6740 0.1880],'HatchDensity',13); 
hatchfill2(object_h(end), 'single','HatchAngle',120,'HatchColor',[0 0.4470 0.7410],'HatchDensity',20);
object_h(end-1).FaceAlpha = 0.15;
object_h(end).FaceAlpha = 0.15;
hMnnc.LineWidth = 1.5;
hMotc.LineWidth = 2;
xlabel('$\kappa$','Interpreter','latex','FontSize',16);
set(gca, 'YGrid', 'on', 'XGrid', 'off')

sgtitle('$\vert {\mathbf{s}_{\mathcal{P}} - \mathbf{s}_{\mathcal{Q}}} \vert_2 = 2 \mathrm{m}$','Interpreter','latex');
set(gcf,'Position',[100 100 900 350]);
saveas(gcf, '/Users/aarongeldert/Documents/MATLAB/THESIS/Figures/ThesisNewFigures/compare2m.png');

%% Cuboid 
load('/Users/aarongeldert/Documents/MATLAB/THESIS/New_ISM_Data/cuboid3_dist2_mapTol1.05/Errors/ERR_cuboid_order3_dist2.mat');

% do same as above; POT works perfectly
% maybe we should do some kind of higher order model with decimation?

%% Error change over distance
clear; clc; close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% Load Data
LD1 = load('/Users/aarongeldert/Documents/MATLAB/THESIS/New_ISM_Data/p1_trapezoidal3_dist0.25_mapTol1.05/Errors/ERR_trapezoidal_order3_dist0.25.mat');
LD2 = load('/Users/aarongeldert/Documents/MATLAB/THESIS/New_ISM_Data/p2_trapezoidal3_dist0.25_mapTol1.05/Errors/ERR_trapezoidal_order3_dist0.25.mat');
LD3 = load('/Users/aarongeldert/Documents/MATLAB/THESIS/New_ISM_Data/trapezoidal3_dist0.5_mapTol1.05/Errors/ERR_trapezoidal_order3_dist0.5.mat');
LD4 = load('/Users/aarongeldert/Documents/MATLAB/THESIS/New_ISM_Data/trapezoidal3_dist1_mapTol1.05/Errors/ERR_trapezoidal_order3_dist1.mat');
LD5 = load('/Users/aarongeldert/Documents/MATLAB/THESIS/New_ISM_Data/trapezoidal3_dist2_mapTol1.05/Errors/ERR_trapezoidal_order3_dist2.mat');
LD6 = load('/Users/aarongeldert/Documents/MATLAB/THESIS/New_ISM_Data/p1_trapezoidal3_dist4_mapTol1.05/Errors/ERR_trapezoidal_order3_dist4.mat');

dists = [0.25 0.5 1 2 4]; % distances
kInd = find(LD1.kVec == 0.5); % when is k = 0.5?

%% get the median error at k=0.5 for each method and distance
close all;
prcTile = 50;
med025 = [prctile([LD0a.ERR.lic(kInd,:),LD0b.ERR.lic(kInd,:)],prcTile,2);...
        prctile([LD0a.ERR.alc(kInd,:),LD0b.ERR.alc(kInd,:)],prcTile,2);...
        prctile([LD0a.ERR.nnc(kInd,:),LD0b.ERR.nnc(kInd,:)],prcTile,2);...
        prctile([LD0a.ERR.otc(kInd,:),LD0b.ERR.otc(kInd,:)],prcTile,2)];
med05 = [prctile(LD0.ERR.lic(kInd,:),prcTile,2);...
        prctile(LD0.ERR.alc(kInd,:),prcTile,2);...
        prctile(LD0.ERR.nnc(kInd,:),prcTile,2);...
        prctile(LD0.ERR.otc(kInd,:),prcTile,2)];
med1 = [prctile(LD1.ERR.lic(kInd,:),prcTile,2);...
        prctile(LD1.ERR.alc(kInd,:),prcTile,2);...
        prctile(LD1.ERR.nnc(kInd,:),prcTile,2);...
        prctile(LD1.ERR.otc(kInd,:),prcTile,2)];
med2 = [prctile(LD2.ERR.lic(kInd,:),prcTile,2);...
        prctile(LD2.ERR.alc(kInd,:),prcTile,2);...
        prctile(LD2.ERR.nnc(kInd,:),prcTile,2);...
        prctile(LD2.ERR.otc(kInd,:),prcTile,2)];
med4 = [prctile([LD4a.ERR.lic(kInd,:),LD4b.ERR.lic(kInd,:)],prcTile,2);...
        prctile([LD4a.ERR.alc(kInd,:),LD4b.ERR.alc(kInd,:)],prcTile,2);...
        prctile([LD4a.ERR.nnc(kInd,:),LD4b.ERR.nnc(kInd,:)],prcTile,2);...
        prctile([LD4a.ERR.otc(kInd,:),LD4b.ERR.otc(kInd,:)],prcTile,2)];

plotData = [med025, med05, med1, med2, med4];

figure(2);
% subplot(121);
loglog(dists, plotData(1,:), '-d','Color',[0.6350 0.0780 0.1840]); hold on;
plot(dists, plotData(2,:), '-*','Color',[0.8500 0.3250 0.098]);
plot(dists, plotData(3,:), '-s','Color',[0.4660 0.6740 0.1880]);
plot(dists, plotData(4,:), '-o','Color',[0 0.4470 0.7410]);

legend('Linear','Aligned Linear','Nearest Neighbor','Partial OT','Location','northwest');
xlabel('$\mathrm{Distance} \ \vert {\mathbf{s}_{\mathcal{P}} - \mathbf{s}_{\mathcal{Q}}} \vert_2 \ \mathrm{[m]}$','Interpreter','latex','FontSize',13); xticks(dists); xlim([min(dists)*0.9 max(dists)*1.1])
ylabel('$\mathrm{Median \ error} \ \mathcal{E}$','Interpreter','latex','FontSize',13);
grid on;

set(gcf,'Position',[100 100 600 400]);
% saveas(gcf, '/Users/aarongeldert/Documents/MATLAB/THESIS/Figures/ThesisNewFigures/compareDists50.png');
