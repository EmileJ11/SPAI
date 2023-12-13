function figgenFunc_VSmaps(params, fignums)
% Created by: Aaron Geldert
% Last modified: 29 Nov 2022
% Generates figures:
% 3: projection map of ground truth VS movement
% 4: projection map of OT interpolation
% 5: projection map of NN interpolation

ISM = params.ISM;
xi = params.xi;
fs = ISM.fs;
geomFile = params.geomFile;
rcvPos1 = params.rcvPos1;
rcvPos2 = params.rcvPos2;

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
figure(fignums(1));

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

% plot continuous ground truth interpolation
for kInd = 2:(numSmooth-1)
    k = kVecSmooth(kInd);
    rcvPos = (1-k)*rcvPos1 + k*rcvPos2;
    PCk = ISM.imageSourcesForRcvPos(rcvPos);
    PCksph = deg2rad(C2S(PCk.pos));
    PCkhap = hammerAidhofProjection(PCksph(:,1:2));
    scatter(PCkhap(:,1), PCkhap(:,2), minDotSize + 0.1*db(1+massScale*PCk.mass).^1.5, [(0.5+0.5*k) .5 (1-0.5*k)], 'filled', 'MarkerFaceAlpha', 1);
end

% generate endpoints + midpoint
for kInd = 1:numSparse
    k = kVecSparse(kInd);
    rcvPos(kInd,:) = (1-k)*rcvPos1 + k*rcvPos2;
    PCgt(kInd) = ISM.imageSourcesForRcvPos(rcvPos(kInd,:));
end

% make struct of PCs
PC1 = PCgt(1);
PC2 = PCgt(end);
PCk = PCgt(2);

% Plot VSs
PCksph = deg2rad(C2S(PCk.pos));
PCkhap = hammerAidhofProjection(PCksph(:,1:2));
hk = scatter(PCkhap(:,1), PCkhap(:,2), minDotSize + db(1+massScale*PCk.mass).^1.5, [.5 .2 .5], 'filled', 'MarkerFaceAlpha', vsAlpha);

PC1sph = deg2rad(C2S(PC1.pos));
PC1hap = hammerAidhofProjection(PC1sph(:,1:2));
h1 = scatter(PC1hap(:,1), PC1hap(:,2), minDotSize + db(1+massScale*PC1.mass).^1.5, [0 0 1], 'filled', 'MarkerFaceAlpha', vsAlpha-0.1);

PC2sph = deg2rad(C2S(PC2.pos));
PC2hap = hammerAidhofProjection(PC2sph(:,1:2));
h2 = scatter(PC2hap(:,1), PC2hap(:,2), minDotSize + db(1+massScale*PC2.mass).^1.5, [1 0 0], 'filled', 'MarkerFaceAlpha', vsAlpha);

xlabel('Azimuth ($^{\circ}$)'); xticks([]);
ylabel('Elevation ($^{\circ}$)'); yticks([]);
set(gcf, 'Position', [100 100 820 500]);

legend([h1 hk h2],'$\mathcal{P}$','$\mathcal{R} \ (\kappa = 0.5)$','$\mathcal{Q}$', 'Location','southwest');

title(['\textbf{' geomName ' room, ground truth}']);

%% Projection map of OT interpolation

figure(fignums(2));

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

% Define POT mapping
OT = PCOT(PCgt(1), PCgt(end), xi);
% OT = VSOT(PCgt(1), PCgt(end), norm(params.rcvPos1 - params.rcvPos2), 0.05);

% define smooth interpolation param values for connection
for kInd = 2:(numSmooth-1)
    k = kVecSmooth(kInd);
    PCk = OT.interpPC(k);
    PCksph = deg2rad(C2S(PCk.pos));
    PCkhap = hammerAidhofProjection(PCksph(:,1:2));
    scatter(PCkhap(:,1), PCkhap(:,2), minDotSize + 0.1*db(1+massScale*PCk.mass).^1.5, [(0.5+0.5*k) .5 (1-0.5*k)], 'filled', 'MarkerFaceAlpha', 1);
end

% make PCs for endpoints and midpoint
PC1 = ISM.imageSourcesForRcvPos(rcvPos1);
PC2 = ISM.imageSourcesForRcvPos(rcvPos2);
PCk = OT.interpPC(0.5);

% Plot VSs
PCksph = deg2rad(C2S(PCk.pos));
PCkhap = hammerAidhofProjection(PCksph(:,1:2));
hk = scatter(PCkhap(:,1), PCkhap(:,2), minDotSize+db(1+massScale*PCk.mass).^1.5, [.5 .2 .5], 'filled', 'MarkerFaceAlpha', vsAlpha);

PC1sph = deg2rad(C2S(PC1.pos));
PC1hap = hammerAidhofProjection(PC1sph(:,1:2));
h1 = scatter(PC1hap(:,1), PC1hap(:,2), minDotSize+db(1+massScale*PC1.mass).^1.5, [0 0 1], 'filled', 'MarkerFaceAlpha', vsAlpha-0.1);

PC2sph = deg2rad(C2S(PC2.pos));
PC2hap = hammerAidhofProjection(PC2sph(:,1:2));
h2 = scatter(PC2hap(:,1), PC2hap(:,2), minDotSize+db(1+massScale*PC2.mass).^1.5, [1 0 0], 'filled', 'MarkerFaceAlpha', vsAlpha);

xlabel('Azimuth ($^{\circ}$)'); xticks([]);
ylabel('Elevation ($^{\circ}$)'); yticks([]);
set(gcf, 'Position', [100 100 820 500]);

legend([h1 hk h2],'$\mathcal{P}$','$\widehat{\mathcal{R}} \ (\kappa = 0.5)$','$\mathcal{Q}$', 'Location','southwest');

title(['\textbf{' geomName ' room, Partial OT mapping method}']);

%% Projection map of NN interpolation
figure(fignums(3));

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


% define smooth interpolation param values for connection
for kInd = 2:(numSmooth-1)
    k = kVecSmooth(kInd);
    PCk = OT.greedyInterp(k);
    PCksph = deg2rad(C2S(PCk.pos));
    PCkhap = hammerAidhofProjection(PCksph(:,1:2));
    scatter(PCkhap(:,1), PCkhap(:,2), minDotSize + 0.1*db(1+massScale*PCk.mass).^1.5, [(0.5+0.5*k) .5 (1-0.5*k)], 'filled', 'MarkerFaceAlpha', 1);
end

% make PCs for NN midpoint
PCk = OT.greedyInterp(0.5);

% Plot VSs
PCksph = deg2rad(C2S(PCk.pos));
PCkhap = hammerAidhofProjection(PCksph(:,1:2));
hk = scatter(PCkhap(:,1), PCkhap(:,2), minDotSize+db(1+massScale*PCk.mass).^1.5, [.5 .2 .5], 'filled', 'MarkerFaceAlpha', vsAlpha);

PC1sph = deg2rad(C2S(PC1.pos));
PC1hap = hammerAidhofProjection(PC1sph(:,1:2));
h1 = scatter(PC1hap(:,1), PC1hap(:,2), minDotSize+db(1+massScale*PC1.mass).^1.5, [0 0 1], 'filled', 'MarkerFaceAlpha', vsAlpha-0.1);

PC2sph = deg2rad(C2S(PC2.pos));
PC2hap = hammerAidhofProjection(PC2sph(:,1:2));
h2 = scatter(PC2hap(:,1), PC2hap(:,2), minDotSize+db(1+massScale*PC2.mass).^1.5, [1 0 0], 'filled', 'MarkerFaceAlpha', vsAlpha);

xlabel('Azimuth ($^{\circ}$)'); xticks([]);
ylabel('Elevation ($^{\circ}$)'); yticks([]);
set(gcf, 'Position', [100 100 820 500]);

legend([h1 hk h2],'$\mathcal{P}$','$\widehat{\mathcal{R}} \ (\kappa = 0.5)$','$\mathcal{Q}$', 'Location','southwest');

title(['\textbf{' geomName ' room, greedy mapping method}']);

end

