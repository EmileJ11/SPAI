%% Figure generation for interpolated projection map
% Created by: Aaron Geldert
% Last modified: 7 Dec 2022

clear; close all; clc;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

% GA Parameters
geomFile = 'trapezoidal.txt';   % which room geometry to use
maxIsOrder = 3;                 % image method maximum order
rCoeff = 0.707;                 % reflection coeff of all surfaces
interpDist = 2;                 % distance the rcv moves between P and Q
wallOffset = 0.5;               % gap btwn a wall and valid src/rcv position
global globalSurfaceCount

% OT Parameters
mapTol = 1.05;                  % cost (scaled to avoid equivalent dummy mappings)
xi = mapTol * interpDist.^2;    % dummy mapping cost
fs = 48e3;                      % sample rate Hz

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

% find 2 random points inside
rcvPos1 = NaN(size(1,3));
rcvPos2 = NaN(size(1,3));
while ~ISM.isPointInsideGeometry(rcvPos1, wallOffset) || ~ISM.isPointInsideGeometry(rcvPos2, wallOffset)
    rcvPos1 = ISM.furthestPt .* rand(1,3);
%     dir = randomPointOnSphere; % random in 3D
    dir = S2C([360*rand(1), 0, 1]); % random in 2D (horizontal)
    rcvPos2 = rcvPos1 + interpDist*dir;
end

geomName = erase(geomFile,'.txt');
geomName(1) = upper(geomName(1));

aziGrid = [0 30  60 90 120 150 180 -30 -60  -120 -150 -90 179.9]'/180*pi;
eleGrid = [-60 -30 0 30 60 ]'/180*pi;
linegray = 0.75; % grid line grey color
vsAlpha = 0.5; % transparency of dots

% smooth connection param values
numSmooth = 501;
kVecSmooth = linspace(0, 1, numSmooth);
% lineScale = 4;
minDotSize = 0.02;
massScale = 40; 

% define sparse (endpts + midpoint) param values for PC plots
kVecSparse = [0 0.5 1];
numSparse = numel(kVecSparse);
rcvPos = NaN(numSparse,3);
PCgt = struct('n',cell(1,numSparse),...
                'pos',cell(1,numSparse),...
                'mass',cell(1,numSparse));
            

%% Projection map of ground truth VS movement
% figure(1);
% 
% % make projection grid
% for ii=1:length(aziGrid)
%     ele=linspace(-pi/2, pi/2, 50);
%     azi=ones(size(ele))*aziGrid(ii);
%     xyGrid = hammerAidhofProjection([azi', ele']);
%     plot(xyGrid(:, 1), xyGrid(:, 2),'--', 'color', [1 1 1]*linegray);
%     hold on
%     text(xyGrid(25, 1),0, num2str(round(180/pi*aziGrid(ii))), ...
%         'fontsize', 9, 'color', [1 1 1]*linegray );
% end
% for jj=1:length(eleGrid)
%     azi=linspace(-pi, pi-.01, 50);
%     ele=ones(size(azi))*eleGrid(jj);
%     xyGrid = hammerAidhofProjection([azi', ele']);
%     plot(xyGrid(:, 1), xyGrid(:, 2),'--', 'color', [1 1 1]*linegray);
%     hold on
%     text(0, xyGrid(25, 2), num2str(180/pi*eleGrid(jj)), ...
%         'fontsize', 9, 'color', [1 1 1]*linegray );
% end
% 
% xlabel('Azimuth ($^{\circ}$)'); xticks([]);
% ylabel('Elevation ($^{\circ}$)'); yticks([]);
% set(gcf, 'Position', [100 100 820 500]);
% title(['\textbf{' geomName ' room, ground truth}']);
% 
% % generate endpoints + midpoint
% for kInd = 1:numSparse
%     k = kVecSparse(kInd);
%     rcvPos(kInd,:) = (1-k)*rcvPos1 + k*rcvPos2;
%     PCgt(kInd) = ISM.imageSourcesForRcvPos(rcvPos(kInd,:));
% end
% % make struct of PCs
% PC1 = PCgt(1);
% PC2 = PCgt(end);
% PCk = PCgt(2);
% 
% % Prepare .gif stuff
% gifname = 'proj_trapezoidal2_gt_ex1.gif';
% tgif = 1/30;
% 
% % plot one endpoint
% PC1sph = deg2rad(C2S(PC1.pos));
% PC1hap = hammerAidhofProjection(PC1sph(:,1:2));
% h1 = scatter(PC1hap(:,1), PC1hap(:,2), minDotSize + db(1+massScale*PC1.mass).^1.5, [0 0 1], 'filled', 'MarkerFaceAlpha', vsAlpha-0.1);
% 
% frame = getframe(gcf); 
% im = frame2im(frame); 
% [imind,cm] = rgb2ind(im,256); 
% imwrite(imind,cm,gifname,'gif', 'Loopcount',inf,'DelayTime',tgif); 
% 
% % plot continuous ground truth interpolation
% for kInd = 2:(numSmooth-1)
%     k = kVecSmooth(kInd);
%     rcvPos = (1-k)*rcvPos1 + k*rcvPos2;
%     PCk = ISM.imageSourcesForRcvPos(rcvPos);
%     PCksph = deg2rad(C2S(PCk.pos));
%     PCkhap = hammerAidhofProjection(PCksph(:,1:2));
%     scatter(PCkhap(:,1), PCkhap(:,2), minDotSize + 0.1*db(1+massScale*PCk.mass).^1.5, [(0.5+0.5*k) .5 (1-0.5*k)], 'filled', 'MarkerFaceAlpha', 1);
%     htemp = scatter(PCkhap(:,1), PCkhap(:,2), minDotSize + db(1+massScale*PCk.mass).^1.5, [(0.5+0.5*k) .5 (1-0.5*k)], 'filled', 'MarkerFaceAlpha', 1);
%     title(['$\kappa = \:$' num2str(round(k,2))], 'FontSize', 15);
%     
% %     pause(0.01);
%     if mod(kInd,3)==0
%         frame = getframe(gcf); 
%         im = frame2im(frame); 
%         [imind,cm] = rgb2ind(im,256); 
%         imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',tgif); 
%     end
%     
%     delete(htemp);
% end
% 
% % plot other endpoint
% PC2sph = deg2rad(C2S(PC2.pos));
% PC2hap = hammerAidhofProjection(PC2sph(:,1:2));
% h2 = scatter(PC2hap(:,1), PC2hap(:,2), minDotSize + db(1+massScale*PC2.mass).^1.5, [1 0 0], 'filled', 'MarkerFaceAlpha', vsAlpha);
% % 
% % % Plot VSs
% % PCksph = deg2rad(C2S(PCk.pos));
% % PCkhap = hammerAidhofProjection(PCksph(:,1:2));
% % hk = scatter(PCkhap(:,1), PCkhap(:,2), minDotSize + db(1+massScale*PCk.mass).^1.5, [.5 .2 .5], 'filled', 'MarkerFaceAlpha', vsAlpha);
% 
% frame = getframe(gcf); 
% im = frame2im(frame); 
% [imind,cm] = rgb2ind(im,256); 
% imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',tgif); 
% 
% % legend([h1 hk h2],'$\mathcal{P}$','$\mathcal{R} \ (\kappa = 0.5)$','$\mathcal{Q}$', 'Location','southwest');


%% First, POT INTERPOLATION:

figure(1);

subplot(5,1,1:4);

% make projection grid
for ii=1:length(aziGrid)
    ele=linspace(-pi/2, pi/2, 50);
    azi=ones(size(ele))*aziGrid(ii);
    xyGrid = hammerAidhofProjection([azi', ele']);
    plot(xyGrid(:, 1), xyGrid(:, 2),'--', 'color', [1 1 1]*linegray);
    hold on
    text(xyGrid(25, 1),0, num2str(round(180/pi*aziGrid(ii))), ...
        'fontsize', 9, 'color', [1 1 1]*linegray );
end
for jj=1:length(eleGrid)
    azi=linspace(-pi, pi-.01, 50);
    ele=ones(size(azi))*eleGrid(jj);
    xyGrid = hammerAidhofProjection([azi', ele']);
    plot(xyGrid(:, 1), xyGrid(:, 2),'--', 'color', [1 1 1]*linegray);
    hold on
    text(0, xyGrid(25, 2), num2str(180/pi*eleGrid(jj)), ...
        'fontsize', 9, 'color', [1 1 1]*linegray );
end

xlabel('Azimuth ($^{\circ}$)'); xticks([]);
ylabel('Elevation ($^{\circ}$)'); yticks([]);
set(gcf, 'Position', [100 100 900 510]);
title('$\kappa = 0$', 'FontSize', 15);

% generate endpoints + midpoint
for kInd = 1:numSparse
    k = kVecSparse(kInd);
    rcvPos(kInd,:) = (1-k)*rcvPos1 + k*rcvPos2;
    PCgt(kInd) = ISM.imageSourcesForRcvPos(rcvPos(kInd,:));
end
% make struct of PCs
PC1 = PCgt(1);
PC2 = PCgt(end);
OT = PCOT(PC1, PC2, xi); % define partial OT plan

% plot one endpoint
PC1sph = deg2rad(C2S(PC1.pos));
PC1hap = hammerAidhofProjection(PC1sph(:,1:2));
h1 = scatter(PC1hap(:,1), PC1hap(:,2), minDotSize + db(1+massScale*PC1.mass).^1.5, [0 0 1]);
hold on;
% plot other endpoint
PC2sph = deg2rad(C2S(PC2.pos));
PC2hap = hammerAidhofProjection(PC2sph(:,1:2));
h2 = scatter(PC2hap(:,1), PC2hap(:,2), minDotSize + db(1+massScale*PC2.mass).^1.5, [1 0 0]);

subplot(5,1,5);
[rir,~] = monoIRfromPC(PC1, 0.04, 4*fs);
stem(db(rir),'BaseValue',-70,'MarkerSize',4,'Color',[0.6 0.6 1]);
hold on;
[rir,~] = monoIRfromPC(PC2, 0.04, 4*fs);
stem(db(rir),'BaseValue',-70,'MarkerSize',4,'Color',[1 0.6 0.6]);
% [rir,~] = monoIRfromPC(PC1, 0.04, 4*fs);
% hs = stem(db(rir),'k','filled','BaseValue',-70,'MarkerSize',3);
axis([1 length(rir) -45 5]);
yticks([-30 0]);
xticks([]);
ylabel('$h(t)$ [dB]');
xlabel('$t$');

disp('Continue?');
pause;

% plot track
subplot(5,1,1:4);
for kInd = 2:(numSmooth-1)
    k = kVecSmooth(kInd);
    rcvPos = (1-k)*rcvPos1 + k*rcvPos2;
%     PCk = ISM.imageSourcesForRcvPos(rcvPos);
    PCk = OT.interpPC(k);
    PCksph = deg2rad(C2S(PCk.pos));
    PCkhap = hammerAidhofProjection(PCksph(:,1:2));
    scatter(PCkhap(:,1), PCkhap(:,2), minDotSize + 0.1*db(1+massScale*PCk.mass).^1.5, [(0.5+0.5*k) .5 (1-0.5*k)], 'filled', 'MarkerFaceAlpha', 1);
end

% Prepare .gif stuff
gifname = ['proj_' geomName '_pot_ex2.gif'];
tgif = 1/30;

frame = getframe(gcf); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 
imwrite(imind,cm,gifname,'gif', 'Loopcount',inf,'DelayTime',tgif); 


% plot moving VSs
for kInd = 2:3:(numSmooth-1)
    subplot(5,1,1:4);
    hold on;
    k = kVecSmooth(kInd);
    rcvPos = (1-k)*rcvPos1 + k*rcvPos2;
%     PCk = ISM.imageSourcesForRcvPos(rcvPos);
    PCk = OT.interpPC(k);
    PCksph = deg2rad(C2S(PCk.pos));
    PCkhap = hammerAidhofProjection(PCksph(:,1:2));
    htemp = scatter(PCkhap(:,1), PCkhap(:,2), minDotSize + db(1+massScale*PCk.mass).^1.5, [(0.5+0.5*k) .5 (1-0.5*k)], 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeColor',[0 0 0]);
    title(['$\kappa = \:$' num2str(round(k,2))], 'FontSize', 15);
    
    subplot(5,1,5);
    hold on;
    [rir, ~] = monoIRfromPC(PCk, 0.04, 4*fs);
    hs = stem(db(rir),'k','filled','BaseValue',-70,'MarkerSize',4,...
        'MarkerFaceColor',[(0.5+0.5*k) .5 (1-0.5*k)], 'MarkerEdgeColor',[0 0 0]);
    axis([1 length(rir) -45 5]);
    yticks([-30 0]);
    xticks([]);
    ylabel('$h(t)$ [dB]');
    xlabel('$t$');
    
%     pause(0.01);
    frame = getframe(gcf); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',tgif); 
    
    delete(htemp);
    delete(hs);
end
% and go back!
for kInd = (numSmooth-1):-3:2
    subplot(5,1,1:4);
    hold on;
    k = kVecSmooth(kInd);
    rcvPos = (1-k)*rcvPos1 + k*rcvPos2;
%     PCk = ISM.imageSourcesForRcvPos(rcvPos);
    PCk = OT.interpPC(k);
    PCksph = deg2rad(C2S(PCk.pos));
    PCkhap = hammerAidhofProjection(PCksph(:,1:2));
    htemp = scatter(PCkhap(:,1), PCkhap(:,2), minDotSize + db(1+massScale*PCk.mass).^1.5, [(0.5+0.5*k) .5 (1-0.5*k)], 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeColor',[0 0 0]);
    title(['$\kappa = \:$' num2str(round(k,2))], 'FontSize', 15);
    
    subplot(5,1,5);
    hold on;
    [rir, ~] = monoIRfromPC(PCk, 0.04, 4*fs);
    hs = stem(db(rir),'k','filled','BaseValue',-70,'MarkerSize',4,...
        'MarkerFaceColor',[(0.5+0.5*k) .5 (1-0.5*k)], 'MarkerEdgeColor',[0 0 0]);
    axis([1 length(rir) -45 5]);
    yticks([-30 0]);
    xticks([]);
    ylabel('$h(t)$ [dB]');
    xlabel('$t$');
    
%     pause(0.01);
    frame = getframe(gcf); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',tgif); 
    
    delete(htemp);
    delete(hs);
end
disp('POT finished. Starting ground truth:');

%% Now, Ground TRuth
figure(2);

subplot(5,1,1:4);

% make projection grid
for ii=1:length(aziGrid)
    ele=linspace(-pi/2, pi/2, 50);
    azi=ones(size(ele))*aziGrid(ii);
    xyGrid = hammerAidhofProjection([azi', ele']);
    plot(xyGrid(:, 1), xyGrid(:, 2),'--', 'color', [1 1 1]*linegray);
    hold on
    text(xyGrid(25, 1),0, num2str(round(180/pi*aziGrid(ii))), ...
        'fontsize', 9, 'color', [1 1 1]*linegray );
end
for jj=1:length(eleGrid)
    azi=linspace(-pi, pi-.01, 50);
    ele=ones(size(azi))*eleGrid(jj);
    xyGrid = hammerAidhofProjection([azi', ele']);
    plot(xyGrid(:, 1), xyGrid(:, 2),'--', 'color', [1 1 1]*linegray);
    hold on
    text(0, xyGrid(25, 2), num2str(180/pi*eleGrid(jj)), ...
        'fontsize', 9, 'color', [1 1 1]*linegray );
end

xlabel('Azimuth ($^{\circ}$)'); xticks([]);
ylabel('Elevation ($^{\circ}$)'); yticks([]);
set(gcf, 'Position', [100 100 900 510]);
title('$\kappa = 0$', 'FontSize', 15);

% generate endpoints + midpoint
for kInd = 1:numSparse
    k = kVecSparse(kInd);
    rcvPos(kInd,:) = (1-k)*rcvPos1 + k*rcvPos2;
    PCgt(kInd) = ISM.imageSourcesForRcvPos(rcvPos(kInd,:));
end
% make struct of PCs
PC1 = PCgt(1);
PC2 = PCgt(end);

% plot one endpoint
PC1sph = deg2rad(C2S(PC1.pos));
PC1hap = hammerAidhofProjection(PC1sph(:,1:2));
h1 = scatter(PC1hap(:,1), PC1hap(:,2), minDotSize + db(1+massScale*PC1.mass).^1.5, [0 0 1]);
hold on;
% plot other endpoint
PC2sph = deg2rad(C2S(PC2.pos));
PC2hap = hammerAidhofProjection(PC2sph(:,1:2));
h2 = scatter(PC2hap(:,1), PC2hap(:,2), minDotSize + db(1+massScale*PC2.mass).^1.5, [1 0 0]);

subplot(5,1,5);
[rir,~] = monoIRfromPC(PC1, 0.04, 4*fs);
stem(db(rir),'BaseValue',-70,'MarkerSize',4,'Color',[0.6 0.6 1]);
hold on;
[rir,~] = monoIRfromPC(PC2, 0.04, 4*fs);
stem(db(rir),'BaseValue',-70,'MarkerSize',4,'Color',[1 0.6 0.6]);
% [rir,~] = monoIRfromPC(PC1, 0.04, 4*fs);
% hs = stem(db(rir),'k','filled','BaseValue',-70,'MarkerSize',3);
axis([1 length(rir) -45 5]);
yticks([-30 0]);
xticks([]);
ylabel('$h(t)$ [dB]');
xlabel('$t$');

% plot track
subplot(5,1,1:4);
for kInd = 2:(numSmooth-1)
    k = kVecSmooth(kInd);
    rcvPos = (1-k)*rcvPos1 + k*rcvPos2;
    PCk = ISM.imageSourcesForRcvPos(rcvPos);
    PCksph = deg2rad(C2S(PCk.pos));
    PCkhap = hammerAidhofProjection(PCksph(:,1:2));
    scatter(PCkhap(:,1), PCkhap(:,2), minDotSize + 0.1*db(1+massScale*PCk.mass).^1.5, [(0.5+0.5*k) .5 (1-0.5*k)], 'filled', 'MarkerFaceAlpha', 1);
end

% Prepare .gif stuff
gifname = ['proj_' geomName '_gt_ex2.gif'];
tgif = 1/30;

frame = getframe(gcf); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 
imwrite(imind,cm,gifname,'gif', 'Loopcount',inf,'DelayTime',tgif); 


% plot moving VSs
for kInd = 2:3:(numSmooth-1)
    subplot(5,1,1:4);
    hold on;
    k = kVecSmooth(kInd);
    rcvPos = (1-k)*rcvPos1 + k*rcvPos2;
    PCk = ISM.imageSourcesForRcvPos(rcvPos);
    PCksph = deg2rad(C2S(PCk.pos));
    PCkhap = hammerAidhofProjection(PCksph(:,1:2));
    htemp = scatter(PCkhap(:,1), PCkhap(:,2), minDotSize + db(1+massScale*PCk.mass).^1.5, [(0.5+0.5*k) .5 (1-0.5*k)], 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeColor',[0 0 0]);
    title(['$\kappa = \:$' num2str(round(k,2))], 'FontSize', 15);
    
    subplot(5,1,5);
    hold on;
    [rir, ~] = monoIRfromPC(PCk, 0.04, 4*fs);
    hs = stem(db(rir),'k','filled','BaseValue',-70,'MarkerSize',4,...
        'MarkerFaceColor',[(0.5+0.5*k) .5 (1-0.5*k)], 'MarkerEdgeColor',[0 0 0]);
    axis([1 length(rir) -45 5]);
    yticks([-30 0]);
    xticks([]);
    ylabel('$h(t)$ [dB]');
    xlabel('$t$');
    
%     pause(0.01);
    frame = getframe(gcf); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',tgif); 
    
    delete(htemp);
    delete(hs);
end
% and go back!
for kInd = (numSmooth-1):-3:2
    subplot(5,1,1:4);
    hold on;
    k = kVecSmooth(kInd);
    rcvPos = (1-k)*rcvPos1 + k*rcvPos2;
    PCk = ISM.imageSourcesForRcvPos(rcvPos);
    PCksph = deg2rad(C2S(PCk.pos));
    PCkhap = hammerAidhofProjection(PCksph(:,1:2));
    htemp = scatter(PCkhap(:,1), PCkhap(:,2), minDotSize + db(1+massScale*PCk.mass).^1.5, [(0.5+0.5*k) .5 (1-0.5*k)], 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeColor',[0 0 0]);
    title(['$\kappa = \:$' num2str(round(k,2))], 'FontSize', 15);
    
    subplot(5,1,5);
    hold on;
    [rir, ~] = monoIRfromPC(PCk, 0.04, 4*fs);
    hs = stem(db(rir),'k','filled','BaseValue',-70,'MarkerSize',4,...
        'MarkerFaceColor',[(0.5+0.5*k) .5 (1-0.5*k)], 'MarkerEdgeColor',[0 0 0]);
    axis([1 length(rir) -45 5]);
    yticks([-30 0]);
    xticks([]);
    ylabel('$h(t)$ [dB]');
    xlabel('$t$');
    
%     pause(0.01);
    frame = getframe(gcf); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',tgif); 
    
    delete(htemp);
    delete(hs);
end
disp('Ground truth done. Starting greedy:');
%% Now, greedy
figure(3);

subplot(5,1,1:4);

% make projection grid
for ii=1:length(aziGrid)
    ele=linspace(-pi/2, pi/2, 50);
    azi=ones(size(ele))*aziGrid(ii);
    xyGrid = hammerAidhofProjection([azi', ele']);
    plot(xyGrid(:, 1), xyGrid(:, 2),'--', 'color', [1 1 1]*linegray);
    hold on
    text(xyGrid(25, 1),0, num2str(round(180/pi*aziGrid(ii))), ...
        'fontsize', 9, 'color', [1 1 1]*linegray );
end
for jj=1:length(eleGrid)
    azi=linspace(-pi, pi-.01, 50);
    ele=ones(size(azi))*eleGrid(jj);
    xyGrid = hammerAidhofProjection([azi', ele']);
    plot(xyGrid(:, 1), xyGrid(:, 2),'--', 'color', [1 1 1]*linegray);
    hold on
    text(0, xyGrid(25, 2), num2str(180/pi*eleGrid(jj)), ...
        'fontsize', 9, 'color', [1 1 1]*linegray );
end

xlabel('Azimuth ($^{\circ}$)'); xticks([]);
ylabel('Elevation ($^{\circ}$)'); yticks([]);
set(gcf, 'Position', [100 100 900 510]);
title('$\kappa = 0$', 'FontSize', 15);

% generate endpoints + midpoint
% for kInd = 1:numSparse
%     k = kVecSparse(kInd);
%     rcvPos(kInd,:) = (1-k)*rcvPos1 + k*rcvPos2;
%     PCgt(kInd) = ISM.imageSourcesForRcvPos(rcvPos(kInd,:));
% end
% make struct of PCs
PC1 = ISM.imageSourcesForRcvPos(rcvPos1);
PC2 = ISM.imageSourcesForRcvPos(rcvPos2);

% plot one endpoint
PC1sph = deg2rad(C2S(PC1.pos));
PC1hap = hammerAidhofProjection(PC1sph(:,1:2));
h1 = scatter(PC1hap(:,1), PC1hap(:,2), minDotSize + db(1+massScale*PC1.mass).^1.5, [0 0 1]);
hold on;
% plot other endpoint
PC2sph = deg2rad(C2S(PC2.pos));
PC2hap = hammerAidhofProjection(PC2sph(:,1:2));
h2 = scatter(PC2hap(:,1), PC2hap(:,2), minDotSize + db(1+massScale*PC2.mass).^1.5, [1 0 0]);

subplot(5,1,5);
[rir,~] = monoIRfromPC(PC1, 0.04, 4*fs);
stem(db(rir),'BaseValue',-70,'MarkerSize',4,'Color',[0.6 0.6 1]);
hold on;
[rir,~] = monoIRfromPC(PC2, 0.04, 4*fs);
stem(db(rir),'BaseValue',-70,'MarkerSize',4,'Color',[1 0.6 0.6]);
% [rir,~] = monoIRfromPC(PC1, 0.04, 4*fs);
% hs = stem(db(rir),'k','filled','BaseValue',-70,'MarkerSize',3);
axis([1 length(rir) -45 5]);
yticks([-30 0]);
xticks([]);
ylabel('$h(t)$ [dB]');
xlabel('$t$');

% plot track
subplot(5,1,1:4);
for kInd = 2:(numSmooth-1)
    k = kVecSmooth(kInd);
    rcvPos = (1-k)*rcvPos1 + k*rcvPos2;
    PCk = OT.greedyInterp(k);
    PCksph = deg2rad(C2S(PCk.pos));
    PCkhap = hammerAidhofProjection(PCksph(:,1:2));
    scatter(PCkhap(:,1), PCkhap(:,2), minDotSize + 0.1*db(1+massScale*PCk.mass).^1.5, [(0.5+0.5*k) .5 (1-0.5*k)], 'filled', 'MarkerFaceAlpha', 1);
end

% Prepare .gif stuff
gifname = ['proj_' geomName '_greedy_ex2.gif'];
tgif = 1/30;

frame = getframe(gcf); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 
imwrite(imind,cm,gifname,'gif', 'Loopcount',inf,'DelayTime',tgif); 


% plot moving VSs
for kInd = 2:3:(numSmooth-1)
    subplot(5,1,1:4);
    hold on;
    k = kVecSmooth(kInd);
    rcvPos = (1-k)*rcvPos1 + k*rcvPos2;
    PCk = OT.greedyInterp(k);
    PCksph = deg2rad(C2S(PCk.pos));
    PCkhap = hammerAidhofProjection(PCksph(:,1:2));
    htemp = scatter(PCkhap(:,1), PCkhap(:,2), minDotSize + db(1+massScale*PCk.mass).^1.5, [(0.5+0.5*k) .5 (1-0.5*k)], 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeColor',[0 0 0]);
    title(['$\kappa = \:$' num2str(round(k,2))], 'FontSize', 15);
    
    subplot(5,1,5);
    hold on;
    [rir, ~] = monoIRfromPC(PCk, 0.04, 4*fs);
    hs = stem(db(rir),'k','filled','BaseValue',-70,'MarkerSize',4,...
        'MarkerFaceColor',[(0.5+0.5*k) .5 (1-0.5*k)], 'MarkerEdgeColor',[0 0 0]);
    axis([1 length(rir) -45 5]);
    yticks([-30 0]);
    xticks([]);
    ylabel('$h(t)$ [dB]');
    xlabel('$t$');
    
%     pause(0.01);
    frame = getframe(gcf); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',tgif); 
    
    delete(htemp);
    delete(hs);
end
% and go back!
for kInd = (numSmooth-1):-3:2
    subplot(5,1,1:4);
    hold on;
    k = kVecSmooth(kInd);
    rcvPos = (1-k)*rcvPos1 + k*rcvPos2;
    PCk = OT.greedyInterp(k);
    PCksph = deg2rad(C2S(PCk.pos));
    PCkhap = hammerAidhofProjection(PCksph(:,1:2));
    htemp = scatter(PCkhap(:,1), PCkhap(:,2), minDotSize + db(1+massScale*PCk.mass).^1.5, [(0.5+0.5*k) .5 (1-0.5*k)], 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeColor',[0 0 0]);
    title(['$\kappa = \:$' num2str(round(k,2))], 'FontSize', 15);
    
    subplot(5,1,5);
    hold on;
    [rir, ~] = monoIRfromPC(PCk, 0.04, 4*fs);
    hs = stem(db(rir),'k','filled','BaseValue',-70,'MarkerSize',4,...
        'MarkerFaceColor',[(0.5+0.5*k) .5 (1-0.5*k)], 'MarkerEdgeColor',[0 0 0]);
    axis([1 length(rir) -45 5]);
    yticks([-30 0]);
    xticks([]);
    ylabel('$h(t)$ [dB]');
    xlabel('$t$');
    
%     pause(0.01);
    frame = getframe(gcf); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',tgif); 
    
    delete(htemp);
    delete(hs);
end