%% Test dataset 1: simulated cuboid room interpolation
% Author: Aaron Geldert
% Last modified: 5 Oct 2022

clear; clc; close all;

tic;
disp('********** STARTING **********');
pathout = '/Users/aarongeldert/Documents/MATLAB/THESIS/Exported Results/IS_movingRcv_staticSrc/';

roomDims = [7.85 5.35 3.15]; % shoebox room dims from ARTSRAM
roomCen = 0.5*roomDims; % diff between the corner and middle of room
roomVertices = [0,              0; 
                roomDims(1),    0;
                roomDims(1),    roomDims(2);
                0,              roomDims(2)];

% PARAMETERS:
srcPos = [6 3.6 1.6]; % source position (fixed)
interpDist = 0.1; % interval between positions (m)
numRcvPos = 21;
mapTol = 2; % OT map tolerance (in xi calculation)
fs = 48e3;

% RCV position, moving:
rcvPos = ones(numRcvPos, 3) .* [2 4 1.8]; % starting position
rcvDir = [1.5 -1 0]; % direction of rcv movement
rcvDir = interpDist .* rcvDir./vecnorm(rcvDir,2);
rcvPos = rcvPos + (0:(numRcvPos-1)).' * rcvDir;

maxOrder = 12;
PCs = struct('n',cell(1,numRcvPos), 'pos',cell(1,numRcvPos), 'mass',cell(1,numRcvPos));

% error parameters
posError = 0.1;
magError = 0.1;

% plot?
doPlot = 0;

% downselection target
maxDist = 0.08* 343; % 80 ms max time
nTarget = 256;
nRange = 0;

for ii = 1:numRcvPos
    tempPC = imageSourcePointCloud(roomDims, srcPos-roomCen, maxOrder);
    
    % add error
    tempPC.mass = tempPC.mass.* (1+ magError*randn(size(tempPC.mass)));
    tempPC.pos =  tempPC.pos + posError*randn(size(tempPC.pos));
    
    % decimate by max distance (100 ms)
    dists = vecnorm(tempPC.pos,2,2);
    directDist = min(dists);
    validInds = dists <= directDist + maxDist;
    tempPC.mass = tempPC.mass(validInds, :);
    tempPC.pos = tempPC.pos(validInds, :);
    disp(['Decimated from ' num2str(tempPC.n) ' to ' num2str(length(tempPC.mass))]);
    
    % image source decimation by mass
    [~, sortInds] = sort(tempPC.mass, 'descend');
    n = nTarget + floor(2*nRange*randn(1) - nRange); 
    if n<10
        n = 10;
    end
    
    sortInds = sortInds(1:n);
    tempPC.mass = tempPC.mass(sortInds, :);
    tempPC.pos = tempPC.pos(sortInds, :) + roomCen - rcvPos(ii,:); % coordinates based on rcvPos
    tempPC.n = n;
    
    PCs(ii) = tempPC; 
end

%% interpolation between pairs
for pairSep = 2:(numRcvPos-1)

    disp(['-> PAIR SEPARATION: ' num2str(pairSep)]);
    ptsPerPair = pairSep - 1;
    numPairs = numRcvPos - pairSep;
    numPts = ptsPerPair * numPairs;
    pairInds = [1:numPairs; (pairSep+1):numRcvPos].';
    
    potError = NaN(numPairs, ptsPerPair);
    linError = NaN(numPairs, ptsPerPair);

    pairDist = pairSep * interpDist;
    xi = mapTol*pairDist; % OT parameter

    scale = 100;
    offset = 10;
    figure(1); 
    parfor ii = 1:numPairs 
        pc1 = PCs(pairInds(ii, 1));
        pc2 = PCs(pairInds(ii, 2));
        disp(['Mapping pair ' num2str(ii)]);
        OT = PCOT(pc1, pc2, xi);
       
        % plot the 2 PCs, room geometry, source and reciever positions 
        cla;
        relSrcPos = srcPos - rcvPos(ii,:);
    
        doPlot = 0;
        if doPlot ~= 0
            roomPoly = polyshape(roomVertices(:,1)-rcvPos(ii,1),...
                roomVertices(:,2)-rcvPos(ii,2));

            plot(roomPoly,'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.1); hold on;
            scatter3(relSrcPos(:,1), relSrcPos(:,2), relSrcPos(:,3), 'sg', 'filled');
            scatter3(0, 0, 0, '.k'); 
            scatter3(pc1.pos(:,1), pc1.pos(:,2), pc1.pos(:,3), offset+ pc1.mass*scale, 'b', 'filled');
            scatter3(pc2.pos(:,1), pc2.pos(:,2), pc2.pos(:,3), offset+ pc2.mass*scale, 'r', 'filled');
            axis image;
            axis([-16 24 -10 15 -20 20]); 
            xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');

            % Plot the OT mappings too
            [row,col,val] = find(OT.T);
            for kk = 1:length(val)
                line([pc1.pos(row(kk),1), pc2.pos(col(kk),1)], [pc1.pos(row(kk),2), pc2.pos(col(kk),2)], [pc1.pos(row(kk),3), pc2.pos(col(kk),3)], 'Color','black');
            end
        end

        for jj = 1:ptsPerPair

            % 3) INTERPOLATION
            % a) Create OT interpolated PC
            k = jj./pairSep;
            [PCot, ~, ~] = OT.interpPC(k);

            % b) Create Linear combination PC
            PClin = combinePCs(pc1, pc2, k);

            % c) Create NN interpolated PC (?)

            % 4) OBJECTIVE EVALUATION
            % Compute OT cost between interpolated PC and ground truth PC
            xiEval = xi; % same number (why? dunno!)
            pcTrue = PCs(pairInds(ii, 1) + jj);

            disp(['Evaluating pt ' num2str(jj)]);
            Ep = getIRerror(PCot, pcTrue, fs, doPlot);
            El = getIRerror(PClin, pcTrue, fs, doPlot);
            potError(ii,jj) = Ep.omni;
            linError(ii,jj) = El.omni;
            
%             potEval = PCOT(PCot, pcTrue, xiEval);
%             linEval = PCOT(PClin, pcTrue, xiEval);
% 
%             potInterpCosts(ii,jj) = sum(potEval.Cx.*potEval.Tx,'all');
%             linInterpCosts(ii,jj) = sum(linEval.Cx.*linEval.Tx,'all');
        end
    end

    filename = [pathout 'IS_movingRcv_staticSrc_sep' num2str(pairSep) '.mat'];
    save(filename, 'PCs', 'potError', 'linError', 'xi', 'pairSep', 'pairInds','posError','magError');
    disp(['<-- Saved ' filename]);
end
toc;

%% Alternatively, make the results here too
ResultsCuboidSim;
