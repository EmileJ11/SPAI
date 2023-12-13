%% Test the Image Source Model
% Created by: Aaron Geldert & Nils Meyer-Kahlen
% Last modified: 24 Oct 2022
clear; close all; clc; 
set(groot,'defaultAxesFontName','Verdana'); % AS PER ICASSP PAPER!
% Prepare the Simulator
global globalSurfaceCount
globalSurfaceCount = 1;

% uh oh, these things are not as in the icassp :/ 
% % ICASSP VERSION
% simulator = IsmSimulator('trapezoidal.txt');
% simulator.maxIsOrder = 2;
% fs = simulator.fs;
% simulator.sourcePosition = [6 2.5 2];
% simulator.simulateImageSources;
% simulator.imageSourcesForRcvPos([2.6 4.3 2]);

% Testing version:
simulator = IsmSimulator('canted.txt');
simulator.maxIsOrder = 2;
fs = simulator.fs;
simulator.sourcePosition = [6 2.5 2];
simulator.simulateImageSources;
simulator.imageSourcesForRcvPos([2.6 4.3 2]);

%% FIGURE 1: Room Geometry
% YZ view (side) 

figure(1); 
subplot(1,5,1:3)
[hs, hr, hi] = simulator.plotGeometry;
legend([hs, hi, hr], 'Source','VS','Receiver','Location','northwest');

view([20 30]);
% set(gca,'YDir','reverse');
% title('Orthographic View');
% axis([-6 11 -6 16 -4 10]);
% xticks([0:5:10]);
% yticks([-5:5:10]);
% zticks([0:5:10]);
grid off;
% set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
% 
% XY view (top)
subplot(1,5,4:5)
simulator.plotGeometry;
view([0 90]);
% set(gca,'YDir','reverse');
% title('Top View');
% axis([-6 11 -6 16 -4 10]);
% xticks([0:5:10]);
% yticks([-5:5:10]);
% zticks([0:5:10]);
% grid off;
xlabel('X');
ylabel('Y');
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
set(gcf,'Position',[100 100 600 340]);

% saveas(gcf, '/Users/aarongeldert/Documents/MATLAB/THESIS/Figures/ThesisFinalFigures/geometry.png');

%% Impulse response comparison plot
% OUTDATED! figgen_RIRalignment.m for the better pairwise comparison plots
k = 0.5;

movDist = 2;
rcvPos1 = [3.2 2.48 1.8];
PC1 = simulator.imageSourcesForRcvPos(rcvPos1); % Pos 1
rcvPos2 = NaN;
% while ~simulator.isPointInsideGeometry(rcvPos2, 0.2)
%     rcvDir = rand(1,3);
%     rcvPos2 = rcvPos1 + movDist*rcvDir./vecnorm(rcvDir);
% end
rcvDir = [0.792207329559554,0.959492426392903,0.655740699156587];
rcvPos2 = rcvPos1 + movDist*rcvDir./vecnorm(rcvDir);
rcvPosk = rcvPos1 + 0.5*movDist*rcvDir./vecnorm(rcvDir);
PC2 = simulator.imageSourcesForRcvPos(rcvPos2); % Pos 2
PCgt = simulator.imageSourcesForRcvPos(rcvPosk); % Pos GT
OT = PCOT(PC1, PC2, (1.5*movDist)^2);
OT2 = VSOT(PC1, PC2, movDist, 0.05);

PClin = combinePCs(PC1, PC2, k); % Linear
PCnn = OT.greedyInterp(k); % NN
PCpot = OT.interpPC(0.5); % POT

%%
sf = hann(round(0.004 * fs));
[ir(:,1), tvec] = monoIRfromPC(PC1, 0.2, fs);
[ir(:,2), ~] = monoIRfromPC(PC2, 0.2, fs);
[ir(:,3), ~] = monoIRfromPC(PCgt, 0.2, fs);
[ir(:,4), ~] = monoIRfromPC(PCpot, 0.2, fs);
[ir(:,5), ~] = monoIRfromPC(PCnn, 0.2, fs);
[ir(:,6), ~] = monoIRfromPC(PClin, 0.2, fs);

% sir = zeros(size(ir));
% for ii = 1:6
%     sir(:,ii) = conv(ir(:,ii),sf,'same');
% end

plotTime = 0.04; % ms of plot desired
dbGap = 40;
dbBln = 35; linBln = 10.^(dbBln/20);
dbOffset = 3;
dbTickGap = -20;
axisLim1 = [0 1e3*plotTime -dbGap-dbBln 0];
axisLim2 = [0 1e3*plotTime -3*dbGap-dbBln 0];
stemWidth = 1.5;

dbIR = db(ir) + dbOffset;
dbIR(dbIR < -dbBln) = -Inf;

greyColor = [.75 .75 .75];

figure(2); clf; h1 = axes;
set(gcf, 'Position', [100 520 600 140]);
set(h1, 'Units', 'Normalized'); % First change to normalized units.
set(h1, 'OuterPosition', [.02, .1, .99, .85]);

% 1: Position P
stem(1e3*tvec, dbIR(:,1),'filled','k',...
    'BaseValue', -dbBln, 'MarkerSize', 3, 'LineWidth', stemWidth);
hold on;
% plot(1e3*tvec, db(sir(:,1))+5,'k');
axis(axisLim1); yticklabels({}); xticklabels({});
set(h1,'YTick',[])

% 2: Position Q
h2 = axes('position',h1.Position);
set(h2,'Color', 'none');
hold on;
stem(1e3*tvec, dbIR(:,2)-dbGap,'filled','k',...
    'BaseValue', -dbBln-dbGap, 'MarkerSize', 3, 'LineWidth', stemWidth);
% plot(1e3*tvec, db(sir(:,2))-dbGap+5,'k');
yticks([-dbGap+dbTickGap -dbGap dbTickGap 0]); 
yticklabels({num2str(dbTickGap), '0', num2str(dbTickGap), '0'});
axis(axisLim1); xticklabels({});
ylabel('$h(t)$ [dB]','Interpreter','latex','FontSize',14);
hold on;
rectangle('Position',[0 -dbBln-dbGap 1e3*plotTime dbBln+dbGap],'FaceColor',[0 0 0 .05],'EdgeColor','none');
hold on;

% labels done manually here
% annotation('textbox', [0.86, 0.85, 0, 0], 'string', '$h_{\mathcal{P}}$','Interpreter','latex',...
%     'FontSize',15);
% annotation('textbox', [0.86, 0.45, 0, 0], 'string', '$h_{\mathcal{Q}}$','Interpreter','latex',...
%     'FontSize',15);

textPosX = 0.16;
annotation('textbox', [textPosX, 0.85, 0, 0], 'string', 'a)',...
    'FontSize',12);
annotation('textbox', [textPosX, 0.50, 0, 0], 'string', 'b)',...
    'FontSize',12);
% xlabel('$t$ (ms)','Interpreter','latex');
% saveas(gcf, '/Users/aarongeldert/Documents/MATLAB/THESIS/Figures/ThesisFinalFigures/irCompare1.png');


% PLOTS 3, 4, 5, 6:
figure(3); clf; h3 = axes;
set(gcf, 'Position', [100 100 600 310]);
set(h3, 'Units', 'Normalized'); % First change to normalized units.
set(h3, 'OuterPosition', [.02, .1, .99, .85]);

% 3: Position R
stem(1e3*tvec, dbIR(:,3),'filled','k',...
    'BaseValue', -dbBln, 'MarkerSize', 3, 'LineWidth', stemWidth);
hold on;
% plot(1e3*tvec, db(sir(:,1))+5,'k');
axis(axisLim2); yticklabels({}); xticklabels({});
set(h3,'YTick',[])

rectangle('Position',[0 -dbBln 1e3*plotTime dbBln],'FaceColor',[0 0 0 .05],'EdgeColor','none');
hold on;

% 4: Estimate POT
h4 = axes('position',h3.Position);
set(h4,'Color', 'none');
hold on;
stem(1e3*tvec, dbIR(:,3)-dbGap,'Color', greyColor,...
    'BaseValue', -dbBln-dbGap, 'MarkerSize', 3, 'LineWidth', 1);
hold on;
stem(1e3*tvec, dbIR(:,4)-dbGap,'filled','Color','#0072BD',...
    'BaseValue', -dbBln-dbGap, 'MarkerSize', 3, 'LineWidth', stemWidth);
% plot(1e3*tvec, db(sir(:,2))-dbGap+5,'k');
axis(axisLim2); xticklabels({}); yticklabels({}); set(h4,'YTick',[]);

% 5: Estimate NN
h5 = axes('position',h4.Position);
set(h5,'Color', 'none');
hold on;
stem(1e3*tvec, dbIR(:,3)-2*dbGap,'Color', greyColor,...
    'BaseValue', -dbBln-2*dbGap, 'MarkerSize', 3, 'LineWidth', 1);
hold on;
stem(1e3*tvec, dbIR(:,5)-2*dbGap,'filled','Color','#77AC30',...
    'BaseValue', -dbBln-2*dbGap, 'MarkerSize', 3, 'LineWidth', stemWidth);
yticks([-3*dbGap+dbTickGap -3*dbGap -2*dbGap+dbTickGap -2*dbGap -dbGap+dbTickGap -dbGap dbTickGap 0]); 
yticklabels({num2str(dbTickGap), '0', num2str(dbTickGap), '0',...
    num2str(dbTickGap), '0',num2str(dbTickGap), '0'});
axis(axisLim2); xticklabels({}); 
ylabel('$h(t)$ [dB]','Interpreter','latex','FontSize',14);

% 6: Estimate Linear
h6 = axes('position',h4.Position);
set(h6,'Color', 'none');
hold on;
stem(1e3*tvec, dbIR(:,3)-3*dbGap,'Color', greyColor,...
    'BaseValue', -dbBln-3*dbGap, 'MarkerSize', 3, 'LineWidth', 1);
hold on;
stem(1e3*tvec, dbIR(:,6)-3*dbGap,'filled','Color','#A2142F',...
    'BaseValue', -dbBln-3*dbGap, 'MarkerSize', 3, 'LineWidth', stemWidth);
hold on;
% plot(1e3*tvec, db(sir(:,3))-2*dbGap+5,'k');
axis(axisLim2); yticklabels({});
set(h6,'YTick',[]);

% labels done manually here
% title('$h_{\mathcal{R}}$','Interpreter','latex');
% annotation('textbox', [0.8, 0.85, 0, 0], 'string', '$h_{\mathcal{R}}$','Interpreter','latex',...
%     'FontSize',15);
% annotation('textbox', [0.8, 0.75, 0, 0], 'string', '$\hat{h}_{\mathcal{R}}$','Interpreter','latex',...
%     'FontSize',15);
% annotation('textbox', [0.8, 0.5, 0, 0], 'string', '$\hat{h}_{\mathcal{R}}$','Interpreter','latex',...
%     'FontSize',15);
% annotation('textbox', [0.8, 0.33, 0, 0], 'string', '$\hat{h}_{\mathcal{R}}$','Interpreter','latex',...
%     'FontSize',15);

annotation('textbox', [textPosX, 0.86, 0, 0], 'string', 'c)',...
    'FontSize',12);
annotation('textbox', [textPosX, 0.71, 0, 0], 'string', 'd)',...
    'FontSize',12);
annotation('textbox', [textPosX, 0.54, 0, 0], 'string', 'e)',...
    'FontSize',12);
annotation('textbox', [textPosX, 0.36, 0, 0], 'string', 'f)',...
    'FontSize',12);
xlabel('$t$ [ms]','Interpreter','latex','FontSize',14);

% saveas(gcf, '/Users/aarongeldert/Documents/MATLAB/THESIS/Figures/ThesisFinalFigures/irCompare2.png');

% sadly, stackedplot does not work for stems :/

