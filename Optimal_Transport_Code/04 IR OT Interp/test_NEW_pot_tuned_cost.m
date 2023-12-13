% Testing if the new cost method shows any improvement with deterministic
% runs
clear; clc; close all;

% Trapezoidal
load('/Users/aarongeldert/Documents/MATLAB/THESIS/EvenNewer_Newest_ISM_Runs/trapezoidal3_dist2_mapTol1.05/Errors/ERR_trapezoidal_order3_dist2.mat');
ERR.otn = [ERR.otn; zeros(1,50)];

axisLim = [0 1 0 0.4];
legendBuff = [-5 -12];

% Setup figure
% set(groot,'defaultAxesFontName','Verdana');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
figure(1); clf; 
subplot(131);

% NEW OT vs pot in DB
[hMlic, hFlic] = iqrPlot(kVec, ERR.otn, 'lin'); hold on;
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');
Legend = {'New POT','Partial OT'};
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