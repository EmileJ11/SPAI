function statisticalCaseResults(geomFile, maxIsOrder, interpDist, numPaths, mapTol, doSave, doPlot)

%% Generate statistical dataset of testing the interpolation methods!
% Created by: Aaron Geldert
% Last modified: 30 Nov 2022

close all;
% tic;
% Setup the room model parameters
% geomFile = 'cuboid.txt';
% geomFile = 'canted.txt';
% geomFile = 'trapezoidal.txt';
fs = 48e3;
% maxIsOrder = 3;
kVec = 0:0.05:1;
numk = length(kVec);
% numPaths = 100; 

% Tunable hyperparameters
% mapTol = 2;
% interpDist = 1; % total distance between pos1 and pos2, meters
rCoeff = 0.707; % about -3db per reflection
wallOffset = 0.5; % source/receiver pos needs to be away from wall

% Initializations
xi = mapTol .* (interpDist.^2); % squared distance of the interpolation gap
ISM = IsmSimulator(geomFile);
ISM.verbose = 0;
ISM.fs = fs;
ISM.maxIsOrder = maxIsOrder;
ISM.rCoeff = rCoeff; 
[~, filenamein, ~] = fileparts(geomFile);

% Make and add paths
parentPath = '/Users/aarongeldert/Documents/MATLAB/THESIS/EvenNewer_Newest_ISM_Runs/';
pathout = [parentPath filenamein num2str(maxIsOrder) '_dist' num2str(interpDist) '_mapTol' num2str(mapTol) '/'];
subpath1 = [pathout 'PCs/'];
subpath2 = [pathout 'Errors/'];
subpath3 = [pathout 'Figures/'];
mkdir(subpath1); 
mkdir(subpath2);  
mkdir(subpath3);
addpath(genpath(pathout));

if ~exist('globalSurfaceCount','var')
    global globalSurfaceCount
end
globalSurfaceCount = 1;

% % Empirically estimate room volume, Sebine RT60
vol.ptsTotal = 1e4;
vol.ptsInner = 0;
vol.rectDims = ISM.furthestPt;
vol.rectVol = prod(vol.rectDims);
for ii = 1:vol.ptsTotal
    pos1 = vol.rectDims .* rand(1,3);
    if ISM.isPointInsideGeometry(pos1)
        vol.ptsInner = vol.ptsInner + 1;
    end
end
vol.realVol = vol.rectVol * (vol.ptsInner/vol.ptsTotal);
vol.rCoeff = rCoeff;
vol.surfaceArea = ISM.totalSurfaceArea;
vol.rt60 = 0.161 * vol.realVol / (vol.surfaceArea * (1-rCoeff)); % 0.5 sec
disp(vol)

%% Generate all the PCs for each path
PCgtc = struct('n',cell(1,numk),...
                'pos',cell(1,numk),...
                'mass',cell(1,numk));
PCotc = PCgtc; % optimal transport clean
PClic = PCgtc; % linear combination clean
PCalc = PCgtc; % aligned linear combination clean
PCnnc = PCgtc; % nearest neighbor clean
PCotn = PCgtc; % optimal transport NEW COST!

% addition of gaussian noise to end positions
% PCgtn = PCgtc;
% PCotn = PCgtc;
% PClin = PCgtc;
% PCaln = PCgtc;
% PCnnn = PCgtc;
            
% init results structs
ERR.otc = zeros(numk, numPaths);
ERR.otn = zeros(numk, numPaths);
ERR.lic = zeros(numk, numPaths);
ERR.alc = zeros(numk, numPaths);
ERR.nnc = zeros(numk, numPaths);

% ERR.otn = zeros(numk, numPaths);
% ERR.lin = zeros(numk, numPaths);
% ERR.aln = zeros(numk, numPaths);
% ERR.nnn = zeros(numk, numPaths);

for nPath = 1:numPaths
    % each path is 1 srcPos, numk rcvPos, and makes interpolated PCs
    disp(['****** PATH ' num2str(nPath) ' of ' num2str(numPaths) ' ******']);
    
    % reset the counter
    globalSurfaceCount = 1;
    
    % Source position
    srcPos = NaN(size(1,3));
    while ~ISM.isPointInsideGeometry(srcPos, wallOffset)
        srcPos = ISM.furthestPt .* rand(1,3); % find a random point inside
    end
    
    % now we have a valid path of interpolation inside the room geometry!
    ISM.sourcePosition = srcPos;
    ISM.simulateImageSources;
    
    % Receiver positions
    % find 2 random points inside
    pos1 = NaN(size(1,3));
    pos2 = NaN(size(1,3));
    while ~ISM.isPointInsideGeometry(pos1, wallOffset) || ~ISM.isPointInsideGeometry(pos2, wallOffset)
        pos1 = ISM.furthestPt .* rand(1,3);
        dir = randomPointOnSphere;
        pos2 = pos1 + interpDist*dir;
    end
    
    rcvPos = NaN(numk,3);
    
    % define the PCs for all k
    for kInd = 1:numk
        k = kVec(kInd);
        rcvPos(kInd,:) = (1-k)*pos1 + k*pos2;
        % ground truth: clean (exact IS model)
        PCgtc(kInd) = ISM.imageSourcesForRcvPos(rcvPos(kInd,:));
        
        % noisy: addition of positional noise, mass variation, mismatching
        % point decimation
%         PCgtn(kInd) = addNoiseToPC(PCgtc(kInd), posNoise, magNoise);
%         PCgtn(kInd) = decimatePC(PCgtn(kInd), nTarget, 2);
    end
    
    % define OT map from endpoints
    
    OTc = PCOT(PCgtc(1), PCgtc(end), xi);
    OTn = VSOT(PCgtc(1), PCgtc(end), interpDist, 0.05);
%     OTn = PCOT(PCgtn(1), PCgtn(end), xi);

    for kInd = 2:(numk-1)
        % interpolate each intermediate k value
        % yes, the endpoints (k=0, k=1) are not set, trivial
        k = kVec(kInd);
        % linear 
        PClic(kInd) = combinePCs(PCgtc(1), PCgtc(end), k);
        PCalc(kInd) = combineAlignedPCs(PCgtc(1), PCgtc(end), k); % aligned
%         PClin(kInd) = combinePCs(PCgtn(1), PCgtn(end), k);
        
        % OT
        PCotc(kInd) = OTc.interpPC(k);
        PCotn(kInd) = OTn.interpPC(k);
        
        % Nearest Neighbor (to be replaced with OT-based version)
%         PCnnc(kInd) = nearestNeighborInterp(PCgtc(1), PCgtc(end), k);
%         PCnnn(kInd) = nearestNeighborInterp(PCgtn(1), PCgtn(end), k);
        
        % Greedy NN algorithm (Oct 24)
        PCnnc(kInd) = OTc.greedyInterp(k);

        % calculate error
        ERR.otc(kInd, nPath) = getIRerror(PCotc(kInd), PCgtc(kInd), fs);
        ERR.otn(kInd, nPath) = getIRerror(PCotn(kInd), PCgtc(kInd), fs);
        ERR.alc(kInd, nPath) = getIRerror(PCalc(kInd), PCgtc(kInd), fs); % aligned
        ERR.lic(kInd, nPath) = getIRerror(PClic(kInd), PCgtc(kInd), fs);
        ERR.nnc(kInd, nPath) = getIRerror(PCnnc(kInd), PCgtc(kInd), fs);
        
%         ERR.otn(kInd, nPath) = getIRerror(PCotn(kInd), PCgtn(kInd), fs); % noisy ground truth!
%         ERR.lin(kInd, nPath) = getIRerror(PClin(kInd), PCgtn(kInd), fs); % noisy ground truth!
%         ERR.aln(kInd, nPath) = getIRerror(PCaln(kInd), PCgtn(kInd), fs); % noisy ground truth!
%         ERR.nnn(kInd, nPath) = getIRerror(PCnnn(kInd), PCgtn(kInd), fs); % noisy ground truth!
    end
    
    % save the PCs per path
    filenameout = [subpath1 'PC_path' num2str(nPath) '.mat'];
    if doSave ~= 0
        save(filenameout, 'PCgtc', 'PCotc','PClic',...
            'srcPos', 'rcvPos', 'nPath', 'interpDist', 'kVec', 'fs', 'maxIsOrder', 'mapTol');
        disp(['<-- Saved ' filenameout]);
    end
    
end

% save the error data
if doSave ~= 0
    filenameout = [subpath2 'ERR_' filenamein '_order' num2str(maxIsOrder)...
        '_dist' num2str(interpDist) '.mat'];
    save(filenameout, 'ERR', 'interpDist', 'kVec', 'fs', 'maxIsOrder', 'mapTol');
    disp(['<-- Saved ' filenameout]);
end
% toc;
% 
%% plots
if doPlot ~= 0
    figure(1);
    ISM.plotGeometry;
    if doSave ~= 0
        saveas(gcf, [subpath3 'geometryISM.fig']);
    end

    %%
    % load('/Users/aarongeldert/Documents/MATLAB/THESIS/ISM_Data/trapezoidal3_dist2_mapTol1/Errors/ERR_trapezoidal_order3_dist2.mat');
    load('/Users/aarongeldert/Documents/MATLAB/THESIS/New_ISM_Data/trapezoidal3_dist0.5_mapTol1.05/Errors/ERR_trapezoidal_order3_dist0.5.mat');

    % ERR.otc = db(1+ERR.otc);
    % ERR.nnc = db(1+ERR.nnc);
    % ERR.lic = db(1+ERR.lic);

    figure(2);
    clf;
    [hMlic, hFlic] = iqrPlot(kVec, ERR.lic, 'lin'); hold on;
    [hMalc, hFalc] = iqrPlot(kVec, ERR.alc, 'al');
    [hMnnc, hFnnc] = iqrPlot(kVec, ERR.nnc, 'nn');
    [hMotc, hFotc] = iqrPlot(kVec, ERR.otc, 'ot');

    xlabel('$\kappa$','Interpreter','latex','FontSize',16);
    ylabel('$\mathcal{E}$','Interpreter','latex','FontSize',16);
    % title(['Interpolation distance: ' num2str(interpDist) ' m']);

    Legend = {'Linear','Aligned Linear','Nearest Neighbor','Partial OT'};
    [~,object_h,~,~] = legendflex([hFlic, hFalc, hFnnc, hFotc],Legend,...
                            'anchor',{'ne','ne'},...  
                            'xscale',0.55,...
                            'buffer',[-3 -3],... 
                            'fontSize',10);

    hatchfill2(object_h(end-3), 'single','HatchAngle',60,'HatchColor',[0.6350 0.0780 0.1840],'HatchDensity',20);
    hatchfill2(object_h(end-2), 'single','HatchAngle',90,'HatchColor',[0.8500 0.3250 0.0980],'HatchDensity',20); 
    hatchfill2(object_h(end-1), 'single','HatchAngle',0,'HatchColor',[0.4660 0.6740 0.1880],'HatchDensity',13); 
    hatchfill2(object_h(end), 'single','HatchAngle',120,'HatchColor',[0 0.4470 0.7410],'HatchDensity',20);
    object_h(end-3).FaceAlpha = 0.15;
    object_h(end-2).FaceAlpha = 0.15;
    object_h(end-1).FaceAlpha = 0.15;
    object_h(end).FaceAlpha = 0.15;

    % legend([hMlic hMnnc hMotc], 'Linear','Nearest Neighbor','Partial OT');
    ylim([0 0.5]); yticks(0:0.1:0.5);
    xticks([0 0.25 0.5 0.75 1]);
    h = gca;
    % h.YGrid = 'on';

    set(gcf,'Position',[100 100 480 280]);

    hMlic.LineWidth = 1;
    hMalc.LineWidth = 1;
    hMnnc.LineWidth = 1;
    hMotc.LineWidth = 2;
    if doSave ~= 0
        saveas(gcf, [subpath3 'error_comparison.png']);
        saveas(gcf, [subpath3 'error_comparison.fig']);
    end

    % figure(3);
    % cla;
    % [hMlin, hFlin] = iqrPlot(kVec, ERR.lin, 'k-');
    % [hMnnn, hFnnn] = iqrPlot(kVec, ERR.nnn, 'k--');
    % [hMotn, hFotn] = iqrPlot(kVec, ERR.otn, 'k-');
    % 
    % xlabel('k (interpolation parameter)');
    % ylabel('Normalized Error');
    % title('Noisy Image Source Simulation ', ['Interpolation distance: ' num2str(interpDist) ' m, ' geomFile]);
    % hatchfill2(hFlin, 'single','HatchAngle',75,'HatchColor',[.75 .75 .75]);
    % hatchfill2(hFnnn, 'single','HatchAngle',0,'HatchColor',[.75 .75 .75]);
    % hatchfill2(hFotn, 'single','HatchAngle',105,'HatchColor',[.5 .5 .5]);
    % legend([hMlin hMnnn hMotn], 'Linear','Nearest Neighbor','Partial OT');
    % ylim([0 0.5]);
    % 
    % hMlin.LineWidth = 1;
    % hMnnn.LineWidth = 1;
    % hMotn.LineWidth = 2;
    % % saveas(gcf, [subpath3 'error_noisy.png']);
    % saveas(gcf, [subpath3 'error_noisy.fig']);
    % toc;

    end
end
