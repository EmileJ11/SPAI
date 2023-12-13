%% Figure Generation for Objective Results
% Created by: Aaron Geldert
% Last modified: 30 Nov 2022

clear; close all; clc;
parentFolder = '/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs';
addpath(genpath('legendflex'));
addpath(genpath(parentFolder));
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

axisLim = [0 1 0 0.4];
legendBuff = [-5 -5];

%% Figure 1: cuboid, 2 meter
LD = load([parentFolder '/cuboid3_dist2_data/Errors/ERR_cuboid_order3_dist2.mat']);
ERR = LD.ERR;
kVec = LD.kVec;

figure(1); clf;

% lin vs pot in DB
t = tiledlayout(1,3); nexttile;
t.TileSpacing = 'compact';
t.Padding = 'compact';
set(gcf,'Position',[100 100 900 350]);

sgtitle('Cuboid room, $\vert {\mathbf{s}_{\mathcal{P}} - \mathbf{s}_{\mathcal{Q}}} \vert_2 = 2 \ [\mathrm{m}]$','Interpreter','latex');

[hMlic, hFlic] = iqrPlot(kVec, ERR.lic, 'lin'); hold on;
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');
Legend = {'Direct linear','Partial OT'};
[~,object_h,~,~] = legendflex([hFlic, hFotc],Legend,...
                            'anchor',{'ne','ne'},...  
                            'xscale',0.55,...
                            'buffer',legendBuff,... 
                            'fontSize',11);
axis(axisLim); yticks(0:0.1:0.5);
ylabel('$\mathcal{E}$','Interpreter','latex','FontSize',15);
xticks([0 0.25 0.5 0.75 1]);
hatchfill2(object_h(end-1), 'single','HatchAngle',75,'HatchColor',[0.6350 0.0780 0.1840],'HatchDensity',20); 
hatchfill2(object_h(end), 'single','HatchAngle',105,'HatchColor',[0 0.4470 0.7410],'HatchDensity',20);
object_h(end-1).FaceAlpha = 0.15;
object_h(end).FaceAlpha = 0.15;
hMlic.LineWidth = 1.5;
hMotc.LineWidth = 2;
xlabel('$\kappa$','Interpreter','latex','FontSize',16);
set(gca, 'YGrid', 'on', 'XGrid', 'off');
                        
% aligned vs pot
nexttile;
[hMalc, hFalc] = iqrPlot(kVec, ERR.alc, 'al'); hold on;
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');
Legend = {'Aligned linear','Partial OT'};
[~,object_h,~,~] = legendflex([hFalc, hFotc],Legend,...
                            'anchor',{'ne','ne'},...  
                            'xscale',0.55,...
                            'buffer',legendBuff,... 
                            'fontSize',11);
axis(axisLim); yticks(0:0.1:0.5);
xticks([0 0.25 0.5 0.75 1]);
object_h(end-1).FaceAlpha = 0.15;
object_h(end).FaceAlpha = 0.15;
hatchfill2(object_h(end-1), 'single','HatchAngle',75,'HatchColor',[0.8500 0.3250 0.098],'HatchDensity',20); 
hatchfill2(object_h(end), 'single','HatchAngle',105,'HatchColor',[0 0.4470 0.7410],'HatchDensity',20);
hMalc.LineWidth = 1.5;
hMotc.LineWidth = 2;
xlabel('$\kappa$','Interpreter','latex','FontSize',16);
set(gca, 'YGrid', 'on', 'XGrid', 'off');

% nn vs pot
nexttile;
[hMnnc, hFnnc] = iqrPlot(kVec, ERR.nnc, 'nn'); hold on;
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');
Legend = {'Greedy mapping','Partial OT'};
[~,object_h,~,~] = legendflex([hFnnc, hFotc],Legend,...
                            'anchor',{'ne','ne'},...  
                            'xscale',0.55,...
                            'buffer',legendBuff,... 
                            'fontSize',11);
axis(axisLim); yticks(0:0.1:0.5);
xticks([0 0.25 0.5 0.75 1]);
hatchfill2(object_h(end-1), 'single','HatchAngle',75,'HatchColor',[0.4660 0.6740 0.1880],'HatchDensity',20); 
hatchfill2(object_h(end), 'single','HatchAngle',105,'HatchColor',[0 0.4470 0.7410],'HatchDensity',20);
object_h(end-1).FaceAlpha = 0.15;
object_h(end).FaceAlpha = 0.15;
hMnnc.LineWidth = 1.5;
hMotc.LineWidth = 2;
xlabel('$\kappa$','Interpreter','latex','FontSize',16);
set(gca, 'YGrid', 'on', 'XGrid', 'off');

saveas(gcf, '/Users/aarongeldert/Documents/MATLAB/THESIS/Figures/ThesisNewFigures/cuboid_2/err_compare.png');

%% Figure 2: canted, 2 meter
LD = load([parentFolder '/canted3_dist2_data/Errors/ERR_canted_order3_dist2.mat']);
ERR = LD.ERR;
kVec = LD.kVec;

figure(2); clf;

% lin vs pot in DB
t = tiledlayout(1,3); nexttile;
t.TileSpacing = 'compact';
t.Padding = 'compact';
set(gcf,'Position',[100 100 900 350]);
sgtitle('Canted room, $\vert {\mathbf{s}_{\mathcal{P}} - \mathbf{s}_{\mathcal{Q}}} \vert_2 = 2 \ [\mathrm{m}]$','Interpreter','latex');

[hMlic, hFlic] = iqrPlot(kVec, ERR.lic, 'lin'); hold on;
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');
Legend = {'Direct linear','Partial OT'};
[~,object_h,~,~] = legendflex([hFlic, hFotc],Legend,...
                            'anchor',{'ne','ne'},...  
                            'xscale',0.55,...
                            'buffer',legendBuff,... 
                            'fontSize',11);
axis(axisLim); yticks(0:0.1:0.5);
ylabel('$\mathcal{E}$','Interpreter','latex','FontSize',15);
xticks([0 0.25 0.5 0.75 1]);
hatchfill2(object_h(end-1), 'single','HatchAngle',75,'HatchColor',[0.6350 0.0780 0.1840],'HatchDensity',20); 
hatchfill2(object_h(end), 'single','HatchAngle',105,'HatchColor',[0 0.4470 0.7410],'HatchDensity',20);
object_h(end-1).FaceAlpha = 0.15;
object_h(end).FaceAlpha = 0.15;
hMlic.LineWidth = 1.5;
hMotc.LineWidth = 2;
xlabel('$\kappa$','Interpreter','latex','FontSize',16);
set(gca, 'YGrid', 'on', 'XGrid', 'off')
                        
% aligned vs pot
nexttile;
[hMalc, hFalc] = iqrPlot(kVec, ERR.alc, 'al'); hold on;
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');
Legend = {'Aligned linear','Partial OT'};
[~,object_h,~,~] = legendflex([hFalc, hFotc],Legend,...
                            'anchor',{'ne','ne'},...  
                            'xscale',0.55,...
                            'buffer',legendBuff,... 
                            'fontSize',11);
axis(axisLim); yticks(0:0.1:0.5);
xticks([0 0.25 0.5 0.75 1]);
object_h(end-1).FaceAlpha = 0.15;
object_h(end).FaceAlpha = 0.15;
hatchfill2(object_h(end-1), 'single','HatchAngle',75,'HatchColor',[0.8500 0.3250 0.098],'HatchDensity',20); 
hatchfill2(object_h(end), 'single','HatchAngle',105,'HatchColor',[0 0.4470 0.7410],'HatchDensity',20);
hMalc.LineWidth = 1.5;
hMotc.LineWidth = 2;
xlabel('$\kappa$','Interpreter','latex','FontSize',16);
set(gca, 'YGrid', 'on', 'XGrid', 'off')

% nn vs pot
nexttile;
[hMnnc, hFnnc] = iqrPlot(kVec, ERR.nnc, 'nn'); hold on;
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');
Legend = {'Greedy mapping','Partial OT'};
[~,object_h,~,~] = legendflex([hFnnc, hFotc],Legend,...
                            'anchor',{'ne','ne'},...  
                            'xscale',0.55,...
                            'buffer',legendBuff,... 
                            'fontSize',11);
axis(axisLim); yticks(0:0.1:0.5);
xticks([0 0.25 0.5 0.75 1]);
hatchfill2(object_h(end-1), 'single','HatchAngle',75,'HatchColor',[0.4660 0.6740 0.1880],'HatchDensity',20); 
hatchfill2(object_h(end), 'single','HatchAngle',105,'HatchColor',[0 0.4470 0.7410],'HatchDensity',20);
object_h(end-1).FaceAlpha = 0.15;
object_h(end).FaceAlpha = 0.15;
hMnnc.LineWidth = 1.5;
hMotc.LineWidth = 2;
xlabel('$\kappa$','Interpreter','latex','FontSize',16);
set(gca, 'YGrid', 'on', 'XGrid', 'off')

saveas(gcf, '/Users/aarongeldert/Documents/MATLAB/THESIS/Figures/ThesisNewFigures/canted_2/err_compare.png');

%% Figure 3: trapezoidal, 2 meter
LD = load([parentFolder '/trapezoidal3_dist2_data/Errors/ERR_trapezoidal_order3_dist2.mat']);
ERR = LD.ERR;
kVec = LD.kVec;

figure(3); clf;

% lin vs pot in DB
t = tiledlayout(1,3); nexttile;
t.TileSpacing = 'compact';
t.Padding = 'compact';
set(gcf,'Position',[100 100 900 350]);
sgtitle('Trapezoidal room, $\vert {\mathbf{s}_{\mathcal{P}} - \mathbf{s}_{\mathcal{Q}}} \vert_2 = 2 \ [\mathrm{m}]$','Interpreter','latex');

[hMlic, hFlic] = iqrPlot(kVec, ERR.lic, 'lin'); hold on;
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');
Legend = {'Direct linear','Partial OT'};
[~,object_h,~,~] = legendflex([hFlic, hFotc],Legend,...
                            'anchor',{'ne','ne'},...  
                            'xscale',0.55,...
                            'buffer',legendBuff,... 
                            'fontSize',11);
axis(axisLim); yticks(0:0.1:0.5);
ylabel('$\mathcal{E}$','Interpreter','latex','FontSize',15);
xticks([0 0.25 0.5 0.75 1]);
hatchfill2(object_h(end-1), 'single','HatchAngle',75,'HatchColor',[0.6350 0.0780 0.1840],'HatchDensity',20); 
hatchfill2(object_h(end), 'single','HatchAngle',105,'HatchColor',[0 0.4470 0.7410],'HatchDensity',20);
object_h(end-1).FaceAlpha = 0.15;
object_h(end).FaceAlpha = 0.15;
hMlic.LineWidth = 1.5;
hMotc.LineWidth = 2;
xlabel('$\kappa$','Interpreter','latex','FontSize',16);
set(gca, 'YGrid', 'on', 'XGrid', 'off')
                        
% aligned vs pot
nexttile;
[hMalc, hFalc] = iqrPlot(kVec, ERR.alc, 'al'); hold on;
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');
Legend = {'Aligned linear','Partial OT'};
[~,object_h,~,~] = legendflex([hFalc, hFotc],Legend,...
                            'anchor',{'ne','ne'},...  
                            'xscale',0.55,...
                            'buffer',legendBuff,... 
                            'fontSize',11);
axis(axisLim); yticks(0:0.1:0.5);
xticks([0 0.25 0.5 0.75 1]);
object_h(end-1).FaceAlpha = 0.15;
object_h(end).FaceAlpha = 0.15;
hatchfill2(object_h(end-1), 'single','HatchAngle',75,'HatchColor',[0.8500 0.3250 0.098],'HatchDensity',20); 
hatchfill2(object_h(end), 'single','HatchAngle',105,'HatchColor',[0 0.4470 0.7410],'HatchDensity',20);
hMalc.LineWidth = 1.5;
hMotc.LineWidth = 2;
xlabel('$\kappa$','Interpreter','latex','FontSize',16);
set(gca, 'YGrid', 'on', 'XGrid', 'off')

% nn vs pot
nexttile;
[hMnnc, hFnnc] = iqrPlot(kVec, ERR.nnc, 'nn'); hold on;
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');
Legend = {'Greedy mapping','Partial OT'};
[~,object_h,~,~] = legendflex([hFnnc, hFotc],Legend,...
                            'anchor',{'ne','ne'},...  
                            'xscale',0.55,...
                            'buffer',legendBuff,... 
                            'fontSize',11);
axis(axisLim); yticks(0:0.1:0.5);
xticks([0 0.25 0.5 0.75 1]);
hatchfill2(object_h(end-1), 'single','HatchAngle',75,'HatchColor',[0.4660 0.6740 0.1880],'HatchDensity',20); 
hatchfill2(object_h(end), 'single','HatchAngle',105,'HatchColor',[0 0.4470 0.7410],'HatchDensity',20);
object_h(end-1).FaceAlpha = 0.15;
object_h(end).FaceAlpha = 0.15;
hMnnc.LineWidth = 1.5;
hMotc.LineWidth = 2;
xlabel('$\kappa$','Interpreter','latex','FontSize',16);
set(gca, 'YGrid', 'on', 'XGrid', 'off')

saveas(gcf, '/Users/aarongeldert/Documents/MATLAB/THESIS/Figures/ThesisNewFigures/trapezoidal_2/err_compare.png');

%% Figure 4: cuboid, 0.5 meter
axisLim = [0 1 0 0.025];

LD = load([parentFolder '/cuboid3_dist0.5_data/Errors/ERR_cuboid_order3_dist0.5.mat']);
ERR = LD.ERR;
kVec = LD.kVec;

figure(4); clf;

% lin vs pot in DB
t = tiledlayout(1,3); nexttile;
t.TileSpacing = 'compact';
t.Padding = 'compact';
set(gcf,'Position',[100 100 900 350]);

sgtitle('Cuboid room, $\vert {\mathbf{s}_{\mathcal{P}} - \mathbf{s}_{\mathcal{Q}}} \vert_2 = 0.5 \ [\mathrm{m}]$','Interpreter','latex');

[hMlic, hFlic] = iqrPlot(kVec, ERR.lic, 'lin'); hold on;
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');
Legend = {'Direct linear','Partial OT'};
[~,object_h,~,~] = legendflex([hFlic, hFotc],Legend,...
                            'anchor',{'ne','ne'},...  
                            'xscale',0.55,...
                            'buffer',legendBuff,... 
                            'fontSize',11);
axis(axisLim); yticks(0:0.005:0.025);
ylabel('$\mathcal{E}$','Interpreter','latex','FontSize',15);
xticks([0 0.25 0.5 0.75 1]);
hatchfill2(object_h(end-1), 'single','HatchAngle',75,'HatchColor',[0.6350 0.0780 0.1840],'HatchDensity',20); 
hatchfill2(object_h(end), 'single','HatchAngle',105,'HatchColor',[0 0.4470 0.7410],'HatchDensity',20);
object_h(end-1).FaceAlpha = 0.15;
object_h(end).FaceAlpha = 0.15;
hMlic.LineWidth = 1.5;
hMotc.LineWidth = 2;
xlabel('$\kappa$','Interpreter','latex','FontSize',16);
set(gca, 'YGrid', 'on', 'XGrid', 'off');
                        
% aligned vs pot
nexttile;
[hMalc, hFalc] = iqrPlot(kVec, ERR.alc, 'al'); hold on;
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');
Legend = {'Aligned linear','Partial OT'};
[~,object_h,~,~] = legendflex([hFalc, hFotc],Legend,...
                            'anchor',{'ne','ne'},...  
                            'xscale',0.55,...
                            'buffer',legendBuff,... 
                            'fontSize',11);
axis(axisLim); yticks(0:0.005:0.025);
xticks([0 0.25 0.5 0.75 1]);
object_h(end-1).FaceAlpha = 0.15;
object_h(end).FaceAlpha = 0.15;
hatchfill2(object_h(end-1), 'single','HatchAngle',75,'HatchColor',[0.8500 0.3250 0.098],'HatchDensity',20); 
hatchfill2(object_h(end), 'single','HatchAngle',105,'HatchColor',[0 0.4470 0.7410],'HatchDensity',20);
hMalc.LineWidth = 1.5;
hMotc.LineWidth = 2;
xlabel('$\kappa$','Interpreter','latex','FontSize',16);
set(gca, 'YGrid', 'on', 'XGrid', 'off');

% nn vs pot
nexttile;
[hMnnc, hFnnc] = iqrPlot(kVec, ERR.nnc, 'nn'); hold on;
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');
Legend = {'Greedy mapping','Partial OT'};
[~,object_h,~,~] = legendflex([hFnnc, hFotc],Legend,...
                            'anchor',{'ne','ne'},...  
                            'xscale',0.55,...
                            'buffer',legendBuff,... 
                            'fontSize',11);
axis(axisLim); yticks(0:0.005:0.025);
xticks([0 0.25 0.5 0.75 1]);
hatchfill2(object_h(end-1), 'single','HatchAngle',75,'HatchColor',[0.4660 0.6740 0.1880],'HatchDensity',20); 
hatchfill2(object_h(end), 'single','HatchAngle',105,'HatchColor',[0 0.4470 0.7410],'HatchDensity',20);
object_h(end-1).FaceAlpha = 0.15;
object_h(end).FaceAlpha = 0.15;
hMnnc.LineWidth = 1.5;
hMotc.LineWidth = 2;
xlabel('$\kappa$','Interpreter','latex','FontSize',16);
set(gca, 'YGrid', 'on', 'XGrid', 'off');

saveas(gcf, '/Users/aarongeldert/Documents/MATLAB/THESIS/Figures/ThesisNewFigures/cuboid_0.5/err_compare.png');

%% Figure 5: canted, 0.5 meter
LD = load([parentFolder '/canted3_dist0.5_data/Errors/ERR_canted_order3_dist0.5.mat']);
ERR = LD.ERR;
kVec = LD.kVec;

figure(5); clf;

% lin vs pot in DB
t = tiledlayout(1,3); nexttile;
t.TileSpacing = 'compact';
t.Padding = 'compact';
set(gcf,'Position',[100 100 900 350]);
sgtitle('Canted room, $\vert {\mathbf{s}_{\mathcal{P}} - \mathbf{s}_{\mathcal{Q}}} \vert_2 = 0.5 \ [\mathrm{m}]$','Interpreter','latex');

[hMlic, hFlic] = iqrPlot(kVec, ERR.lic, 'lin'); hold on;
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');
Legend = {'Direct linear','Partial OT'};
[~,object_h,~,~] = legendflex([hFlic, hFotc],Legend,...
                            'anchor',{'ne','ne'},...  
                            'xscale',0.55,...
                            'buffer',legendBuff,... 
                            'fontSize',11);
axis(axisLim); yticks(0:0.005:0.025);
ylabel('$\mathcal{E}$','Interpreter','latex','FontSize',15);
xticks([0 0.25 0.5 0.75 1]);
hatchfill2(object_h(end-1), 'single','HatchAngle',75,'HatchColor',[0.6350 0.0780 0.1840],'HatchDensity',20); 
hatchfill2(object_h(end), 'single','HatchAngle',105,'HatchColor',[0 0.4470 0.7410],'HatchDensity',20);
object_h(end-1).FaceAlpha = 0.15;
object_h(end).FaceAlpha = 0.15;
hMlic.LineWidth = 1.5;
hMotc.LineWidth = 2;
xlabel('$\kappa$','Interpreter','latex','FontSize',16);
set(gca, 'YGrid', 'on', 'XGrid', 'off')
                        
% aligned vs pot
nexttile;
[hMalc, hFalc] = iqrPlot(kVec, ERR.alc, 'al'); hold on;
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');
Legend = {'Aligned linear','Partial OT'};
[~,object_h,~,~] = legendflex([hFalc, hFotc],Legend,...
                            'anchor',{'ne','ne'},...  
                            'xscale',0.55,...
                            'buffer',legendBuff,... 
                            'fontSize',11);
axis(axisLim); yticks(0:0.005:0.025);
xticks([0 0.25 0.5 0.75 1]);
object_h(end-1).FaceAlpha = 0.15;
object_h(end).FaceAlpha = 0.15;
hatchfill2(object_h(end-1), 'single','HatchAngle',75,'HatchColor',[0.8500 0.3250 0.098],'HatchDensity',20); 
hatchfill2(object_h(end), 'single','HatchAngle',105,'HatchColor',[0 0.4470 0.7410],'HatchDensity',20);
hMalc.LineWidth = 1.5;
hMotc.LineWidth = 2;
xlabel('$\kappa$','Interpreter','latex','FontSize',16);
set(gca, 'YGrid', 'on', 'XGrid', 'off')

% nn vs pot
nexttile;
[hMnnc, hFnnc] = iqrPlot(kVec, ERR.nnc, 'nn'); hold on;
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');
Legend = {'Greedy mapping','Partial OT'};
[~,object_h,~,~] = legendflex([hFnnc, hFotc],Legend,...
                            'anchor',{'ne','ne'},...  
                            'xscale',0.55,...
                            'buffer',legendBuff,... 
                            'fontSize',11);
axis(axisLim); yticks(0:0.005:0.025);
xticks([0 0.25 0.5 0.75 1]);
hatchfill2(object_h(end-1), 'single','HatchAngle',75,'HatchColor',[0.4660 0.6740 0.1880],'HatchDensity',20); 
hatchfill2(object_h(end), 'single','HatchAngle',105,'HatchColor',[0 0.4470 0.7410],'HatchDensity',20);
object_h(end-1).FaceAlpha = 0.15;
object_h(end).FaceAlpha = 0.15;
hMnnc.LineWidth = 1.5;
hMotc.LineWidth = 2;
xlabel('$\kappa$','Interpreter','latex','FontSize',16);
set(gca, 'YGrid', 'on', 'XGrid', 'off')

saveas(gcf, '/Users/aarongeldert/Documents/MATLAB/THESIS/Figures/ThesisNewFigures/canted_0.5/err_compare.png');

%% Figure 6: trapezoidal, 0.5 meter
LD = load([parentFolder '/trapezoidal3_dist0.5_data/Errors/ERR_trapezoidal_order3_dist0.5.mat']);
ERR = LD.ERR;
kVec = LD.kVec;

figure(6); clf;

% lin vs pot in DB
t = tiledlayout(1,3); nexttile;
t.TileSpacing = 'compact';
t.Padding = 'compact';
set(gcf,'Position',[100 100 900 350]);
sgtitle('Trapezoidal room, $\vert {\mathbf{s}_{\mathcal{P}} - \mathbf{s}_{\mathcal{Q}}} \vert_2 = 0.5 \ [\mathrm{m}]$','Interpreter','latex');

[hMlic, hFlic] = iqrPlot(kVec, ERR.lic, 'lin'); hold on;
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');
Legend = {'Direct linear','Partial OT'};
[~,object_h,~,~] = legendflex([hFlic, hFotc],Legend,...
                            'anchor',{'ne','ne'},...  
                            'xscale',0.55,...
                            'buffer',legendBuff,... 
                            'fontSize',11);
axis(axisLim); yticks(0:0.005:0.025);
ylabel('$\mathcal{E}$','Interpreter','latex','FontSize',15);
xticks([0 0.25 0.5 0.75 1]);
hatchfill2(object_h(end-1), 'single','HatchAngle',75,'HatchColor',[0.6350 0.0780 0.1840],'HatchDensity',20); 
hatchfill2(object_h(end), 'single','HatchAngle',105,'HatchColor',[0 0.4470 0.7410],'HatchDensity',20);
object_h(end-1).FaceAlpha = 0.15;
object_h(end).FaceAlpha = 0.15;
hMlic.LineWidth = 1.5;
hMotc.LineWidth = 2;
xlabel('$\kappa$','Interpreter','latex','FontSize',16);
set(gca, 'YGrid', 'on', 'XGrid', 'off')
                        
% aligned vs pot
nexttile;
[hMalc, hFalc] = iqrPlot(kVec, ERR.alc, 'al'); hold on;
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');
Legend = {'Aligned linear','Partial OT'};
[~,object_h,~,~] = legendflex([hFalc, hFotc],Legend,...
                            'anchor',{'ne','ne'},...  
                            'xscale',0.55,...
                            'buffer',legendBuff,... 
                            'fontSize',11);
axis(axisLim); yticks(0:0.005:0.025);
xticks([0 0.25 0.5 0.75 1]);
object_h(end-1).FaceAlpha = 0.15;
object_h(end).FaceAlpha = 0.15;
hatchfill2(object_h(end-1), 'single','HatchAngle',75,'HatchColor',[0.8500 0.3250 0.098],'HatchDensity',20); 
hatchfill2(object_h(end), 'single','HatchAngle',105,'HatchColor',[0 0.4470 0.7410],'HatchDensity',20);
hMalc.LineWidth = 1.5;
hMotc.LineWidth = 2;
xlabel('$\kappa$','Interpreter','latex','FontSize',16);
set(gca, 'YGrid', 'on', 'XGrid', 'off')

% nn vs pot
nexttile;
[hMnnc, hFnnc] = iqrPlot(kVec, ERR.nnc, 'nn'); hold on;
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');
Legend = {'Greedy mapping','Partial OT'};
[~,object_h,~,~] = legendflex([hFnnc, hFotc],Legend,...
                            'anchor',{'ne','ne'},...  
                            'xscale',0.55,...
                            'buffer',legendBuff,... 
                            'fontSize',11);
axis(axisLim); yticks(0:0.005:0.025);
xticks([0 0.25 0.5 0.75 1]);
hatchfill2(object_h(end-1), 'single','HatchAngle',75,'HatchColor',[0.4660 0.6740 0.1880],'HatchDensity',20); 
hatchfill2(object_h(end), 'single','HatchAngle',105,'HatchColor',[0 0.4470 0.7410],'HatchDensity',20);
object_h(end-1).FaceAlpha = 0.15;
object_h(end).FaceAlpha = 0.15;
hMnnc.LineWidth = 1.5;
hMotc.LineWidth = 2;
xlabel('$\kappa$','Interpreter','latex','FontSize',16);
set(gca, 'YGrid', 'on', 'XGrid', 'off')

saveas(gcf, '/Users/aarongeldert/Documents/MATLAB/THESIS/Figures/ThesisNewFigures/trapezoidal_0.5/err_compare.png');

%% Figure 7: Cuboid room, median error over distance
% Load Data
clear;
LD(1) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/cuboid3_dist0.125_data/Errors/ERR_cuboid_order3_dist0.125.mat');
LD(2) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/cuboid3_dist0.25_data/Errors/ERR_cuboid_order3_dist0.25.mat');
LD(3) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/cuboid3_dist0.5_data/Errors/ERR_cuboid_order3_dist0.5.mat');
LD(4) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/cuboid3_dist1_data/Errors/ERR_cuboid_order3_dist1.mat');
LD(5) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/cuboid3_dist2_data/Errors/ERR_cuboid_order3_dist2.mat');
LD(6) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/cuboid3_dist4_data/Errors/ERR_cuboid_order3_dist4.mat');

dists = [0.125 0.25 0.5 1 2 4];
numDists = numel(dists);
kInd = find(LD(1).kVec == 0.5); % where is midpoint data?
pTile = 50;
meds = NaN(numDists, 4);

for ii = 1:numDists
    meds(ii,1) = prctile(LD(ii).ERR.lic(kInd,:), 50, 2);
    meds(ii,2) = prctile(LD(ii).ERR.alc(kInd,:), 50, 2);
    meds(ii,3) = prctile(LD(ii).ERR.nnc(kInd,:), 50, 2);
    meds(ii,4) = prctile(LD(ii).ERR.otc(kInd,:), 50, 2);
end

figure(7); clf;
loglog(dists, meds(:,1), '-d','Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'MarkerSize',8); hold on;
plot(dists, meds(:,2), '-*','Color',[0.8500 0.3250 0.098],'LineWidth',1.5,'MarkerSize',8);
plot(dists, meds(:,3), '-s','Color',[0.4660 0.6740 0.1880],'LineWidth',1.5,'MarkerSize',8);
plot(dists, meds(:,4), '-o','Color',[0 0.4470 0.7410],'LineWidth',1.5,'MarkerSize',8);

legend('Direct linear','Aligned linear','Greedy mapping','Partial OT','Location','southeast');
xlabel('$\mathrm{Distance} \ \vert {\mathbf{s}_{\mathcal{P}} - \mathbf{s}_{\mathcal{Q}}} \vert_2 \ \mathrm{[m]}$','Interpreter','latex','FontSize',13); 
xticks(dists); xlim([min(dists)*0.9 max(dists)*1.1])
ylim([0.75e-5 1]);
ylabel('$\mathrm{Median \ error} \ \mathcal{E}$','Interpreter','latex','FontSize',13);
grid on; 
grid minor; grid minor;

set(gcf,'Position',[100 100 700 360]);
title('\textbf{Cuboid room}, $\kappa = 0.5$','FontSize',13);
% saveas(gcf, '/Users/aarongeldert/Documents/MATLAB/THESIS/Figures/ThesisNewFigures/cuboid_dist_error.png');


%% Figure 8: Canted room, median error over distance
% Load Data
clear;
LD(1) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/canted3_dist0.125_data/Errors/ERR_canted_order3_dist0.125.mat');
LD(2) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/canted3_dist0.25_data/Errors/ERR_canted_order3_dist0.25.mat');
LD(3) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/canted3_dist0.5_data/Errors/ERR_canted_order3_dist0.5.mat');
LD(4) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/canted3_dist1_data/Errors/ERR_canted_order3_dist1.mat');
LD(5) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/canted3_dist2_data/Errors/ERR_canted_order3_dist2.mat');
LD(6) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/canted3_dist4_data/Errors/ERR_canted_order3_dist4.mat');

dists = [0.125 0.25 0.5 1 2 4];
numDists = numel(dists);
kInd = find(LD(1).kVec == 0.5); % where is midpoint data?
pTile = 50;
meds = NaN(numDists, 4);

for ii = 1:numDists
    meds(ii,1) = prctile(LD(ii).ERR.lic(kInd,:), 50, 2);
    meds(ii,2) = prctile(LD(ii).ERR.alc(kInd,:), 50, 2);
    meds(ii,3) = prctile(LD(ii).ERR.nnc(kInd,:), 50, 2);
    meds(ii,4) = prctile(LD(ii).ERR.otc(kInd,:), 50, 2);
end

figure(8); clf;
loglog(dists, meds(:,1), '-d','Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'MarkerSize',8); hold on;
plot(dists, meds(:,2), '-*','Color',[0.8500 0.3250 0.098],'LineWidth',1.5,'MarkerSize',8);
plot(dists, meds(:,3), '-s','Color',[0.4660 0.6740 0.1880],'LineWidth',1.5,'MarkerSize',8);
plot(dists, meds(:,4), '-o','Color',[0 0.4470 0.7410],'LineWidth',1.5,'MarkerSize',8);

legend('Direct linear','Aligned linear','Greedy mapping','Partial OT','Location','southeast');
xlabel('$\mathrm{Distance} \ \vert {\mathbf{s}_{\mathcal{P}} - \mathbf{s}_{\mathcal{Q}}} \vert_2 \ \mathrm{[m]}$','Interpreter','latex','FontSize',13); 
xticks(dists); xlim([min(dists)*0.9 max(dists)*1.1])
ylim([0.75e-5 1]);
ylabel('$\mathrm{Median \ error} \ \mathcal{E}$','Interpreter','latex','FontSize',13);
grid on; 
grid minor; grid minor;

set(gcf,'Position',[100 100 700 360]);
title('\textbf{Canted room}, $\kappa = 0.5$','FontSize',13);
saveas(gcf, '/Users/aarongeldert/Documents/MATLAB/THESIS/Figures/ThesisNewFigures/canted_dist_error.png');

%% Figure 9: Trapezoidal room, median error over distance
% Load Data
clear;
LD(1) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/trapezoidal3_dist0.125_data/Errors/ERR_trapezoidal_order3_dist0.125.mat');
LD(2) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/trapezoidal3_dist0.25_data/Errors/ERR_trapezoidal_order3_dist0.25.mat');
LD(3) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/trapezoidal3_dist0.5_data/Errors/ERR_trapezoidal_order3_dist0.5.mat');
LD(4) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/trapezoidal3_dist1_data/Errors/ERR_trapezoidal_order3_dist1.mat');
LD(5) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/trapezoidal3_dist2_data/Errors/ERR_trapezoidal_order3_dist2.mat');
LD(6) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/trapezoidal3_dist4_data/Errors/ERR_trapezoidal_order3_dist4.mat');

dists = [0.125 0.25 0.5 1 2 4];
numDists = numel(dists);
kInd = find(LD(1).kVec == 0.5); % where is midpoint data?
pTile = 50;
meds = NaN(numDists, 4);

for ii = 1:numDists
    meds(ii,1) = prctile(LD(ii).ERR.lic(kInd,:), 50, 2);
    meds(ii,2) = prctile(LD(ii).ERR.alc(kInd,:), 50, 2);
    meds(ii,3) = prctile(LD(ii).ERR.nnc(kInd,:), 50, 2);
    meds(ii,4) = prctile(LD(ii).ERR.otc(kInd,:), 50, 2);
end

figure(9); clf;
loglog(dists, meds(:,1), '-d','Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'MarkerSize',8); hold on;
plot(dists, meds(:,2), '-*','Color',[0.8500 0.3250 0.098],'LineWidth',1.5,'MarkerSize',8);
plot(dists, meds(:,3), '-s','Color',[0.4660 0.6740 0.1880],'LineWidth',1.5,'MarkerSize',8);
plot(dists, meds(:,4), '-o','Color',[0 0.4470 0.7410],'LineWidth',1.5,'MarkerSize',8);

legend('Direct linear','Aligned linear','Greedy mapping','Partial OT','Location','southeast');
xlabel('$\mathrm{Distance} \ \vert {\mathbf{s}_{\mathcal{P}} - \mathbf{s}_{\mathcal{Q}}} \vert_2 \ \mathrm{[m]}$','Interpreter','latex','FontSize',13); 
xticks(dists); xlim([min(dists)*0.9 max(dists)*1.1])
ylim([0.75e-5 1]);
ylabel('$\mathrm{Median \ error} \ \mathcal{E}$','Interpreter','latex','FontSize',13);
grid on; 
grid minor; grid minor;

set(gcf,'Position',[100 100 700 360]);
title('\textbf{Trapezoidal room, $\kappa = 0.5$}','FontSize',13);
saveas(gcf, '/Users/aarongeldert/Documents/MATLAB/THESIS/Figures/ThesisNewFigures/trapezoidal_dist_error.png');

%% Ranking plots
clear;
figure();
spi = 1;

LD(1,1) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/cuboid3_dist0.125_data/Errors/ERR_cuboid_order3_dist0.125.mat');
LD(1,2) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/cuboid3_dist0.25_data/Errors/ERR_cuboid_order3_dist0.25.mat');
LD(1,3) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/cuboid3_dist0.5_data/Errors/ERR_cuboid_order3_dist0.5.mat');
LD(1,4) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/cuboid3_dist1_data/Errors/ERR_cuboid_order3_dist1.mat');
LD(1,5) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/cuboid3_dist2_data/Errors/ERR_cuboid_order3_dist2.mat');
LD(1,6) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/cuboid3_dist4_data/Errors/ERR_cuboid_order3_dist4.mat');
LD(2,1) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/canted3_dist0.125_data/Errors/ERR_canted_order3_dist0.125.mat');
LD(2,2) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/canted3_dist0.25_data/Errors/ERR_canted_order3_dist0.25.mat');
LD(2,3) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/canted3_dist0.5_data/Errors/ERR_canted_order3_dist0.5.mat');
LD(2,4) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/canted3_dist1_data/Errors/ERR_canted_order3_dist1.mat');
LD(2,5) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/canted3_dist2_data/Errors/ERR_canted_order3_dist2.mat');
LD(2,6) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/canted3_dist4_data/Errors/ERR_canted_order3_dist4.mat');
LD(3,1) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/trapezoidal3_dist0.125_data/Errors/ERR_trapezoidal_order3_dist0.125.mat');
LD(3,2) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/trapezoidal3_dist0.25_data/Errors/ERR_trapezoidal_order3_dist0.25.mat');
LD(3,3) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/trapezoidal3_dist0.5_data/Errors/ERR_trapezoidal_order3_dist0.5.mat');
LD(3,4) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/trapezoidal3_dist1_data/Errors/ERR_trapezoidal_order3_dist1.mat');
LD(3,5) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/trapezoidal3_dist2_data/Errors/ERR_trapezoidal_order3_dist2.mat');
LD(3,6) = load('/Users/aarongeldert/Documents/MATLAB/THESIS/Newest_ISM_Runs/trapezoidal3_dist4_data/Errors/ERR_trapezoidal_order3_dist4.mat');

for roomInd = 1:3
    distanceData = zeros(4,6);
    for distInd = 1:6
    
        countData = determineRankings(LD(roomInd,distInd).ERR);
        countData = 100*countData./sum(countData); % convert to percentage
        % first column is the best result data
        distanceData(:,distInd) = countData(:,1);
    end
    
    subplot(1,3,spi); spi = spi+1;
    b = bar(1:6, distanceData.', 'stacked');
    b(1).FaceColor = [0 0.4470 0.7410]; % POT
    b(2).FaceColor = [0.4660 0.6740 0.1880]; % Greedy
    b(3).FaceColor = [0.8500 0.3250 0.0980]; % aligned
    b(4).FaceColor = [0.6350 0.0780 0.1840]; % linear
    xticks(1:6);
    xticklabels({'0.125','0.25','0.5','1','2','4'});
    xtickangle(45);
    ylim([0 100]);
    yticks([0:10:100]);
end
subplot(131);
title('Cuboid');
ylabel('Occurrence (\%)');

subplot(132);
title('Canted');
xlabel('Distance (m)');

subplot(133);
title('Trapezoidal');
legend('Partial OT','Greedy mapping','Simple linear','Aligned linear','Location','southeast');

sgtitle('Method with lowest RIR error');
set(gcf,'Position',[100 100 700 300]);
saveas(gcf,'/Users/aarongeldert/Documents/MATLAB/THESIS/Figures/ThesisNewFigures/bestMethodBarPlot.png');