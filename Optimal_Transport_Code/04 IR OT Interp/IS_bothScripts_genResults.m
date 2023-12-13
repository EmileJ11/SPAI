%% BOTH SCRIPTS to produce results
% combines IS_movingSrc_staticRcv.m and IS_movingRcv_staticSrc.m
% -> fixed the problem with interpDist, I hope
% -> results should be reasonable and describe the OT benefits?

% Author: Aaron Geldert
% Last modified: 23 Aug 2022

%% MOVING SRC, STATIC RCV

clear; clc; close all;

tic;
disp('******** MOVING SRC, STATIC RCV ********');
pathout = '/Users/aarongeldert/Documents/MATLAB/THESIS/Exported Results/IS_movingSrc_staticRcv/';

roomDims = [7.85 5.35 3.15]; % shoebox room dims from ARTSRAM
roomCen = 0.5*roomDims;
roomVertices = [0,              0; 
                roomDims(1),    0;
                roomDims(1),    roomDims(2);
                0,              roomDims(2)];

rcvPos = [2 3.3 1.6]; % arbitrary pos, origin at lower left ground corner

interpDist = 0.1; % interval between positions
numSrcPos = 26;

srcPos = ones(numSrcPos, 3) .* [6.2 3.3 1.8]; % starting position
srcDir = [-3 -2 0]; % direction of src movement
srcDir = interpDist .* srcDir./vecnorm(srcDir,2);
srcPos = srcPos + (0:(numSrcPos-1)).' * srcDir;

maxOrder = 10;
PCs = struct('n',cell(1,numSrcPos), 'pos',cell(1,numSrcPos), 'mass',cell(1,numSrcPos));

% error parameters
posError = 0.2;
magError = 0.2;

% downselection target
nTarget = 64;
nRange = 3;

% OT map tolerance (in xi calculation)
mapTol = 4;

for ii = 1:numSrcPos
    
    % note, the geometry of IS keeps changing since src is now moving!
    tempPC = imageSourcePointCloud(roomDims, srcPos(ii,:)-roomCen, maxOrder);
    
    % add error
    tempPC.mass = tempPC.mass .* ((1-magError/2)+ magError*rand(size(tempPC.mass)));
    tempPC.pos =  tempPC.pos .* ((1-posError/2)+ posError*rand(size(tempPC.pos)));

    % image source decimation by mass
    [~, sortInds] = sort(tempPC.mass, 'descend');
    n = nTarget + floor(2*nRange*randn(1) - nRange); 
    if n<10
        n = 10;
    end
    
    sortInds = sortInds(1:n);
    tempPC.mass = tempPC.mass(sortInds, :);
    tempPC.pos = tempPC.pos(sortInds, :) + roomCen - rcvPos; % coordinates based on rcvPos
    tempPC.n = n;
    
    PCs(ii) = tempPC; 
end

% interpolation between pairs
for pairSep = 2:25
    toc;
    disp(['-> PAIR SEPARATION: ' num2str(pairSep)]);
    ptsPerPair = pairSep - 1;
    numPairs = numSrcPos - pairSep;
    numPts = ptsPerPair * numPairs;
    pairInds = [1:numPairs; (pairSep+1):numSrcPos].';
    potInterpCosts = NaN(numPairs, ptsPerPair);
    linInterpCosts = NaN(numPairs, ptsPerPair); 
    
    pairDist = pairSep * interpDist;
    xi = mapTol*pairDist; % OT parameter

    scale = 100;
    offset = 10;
    figure(1); 
    for ii = 1:numPairs 
        pc1 = PCs(pairInds(ii, 1));
        pc2 = PCs(pairInds(ii, 2));
        disp(['Mapping pair ' num2str(ii)]);
        OT = PCOT(pc1, pc2, xi);

        % plot the 2 PCs, room geometry, source and reciever positions 
        cla;
        relSrcPos = srcPos(ii,:) - rcvPos;

        roomPoly = polyshape(roomVertices(:,1)-rcvPos(:,1),...
            roomVertices(:,2)-rcvPos(:,2));

        plot(roomPoly,'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.1); hold on;
        scatter3(relSrcPos(:,1), relSrcPos(:,2), relSrcPos(:,3), '.g');
        scatter3(0, 0, 0, 'sk','filled'); 
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

        parfor jj = 1:ptsPerPair

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
            potEval = PCOT(PCot, pcTrue, xiEval);
            linEval = PCOT(PClin, pcTrue, xiEval);

            potInterpCosts(ii,jj) = sum(potEval.Cx.*potEval.Tx,'all');
            linInterpCosts(ii,jj) = sum(linEval.Cx.*linEval.Tx,'all');
        end
    end

    filename = [pathout 'IS_movingSrc_staticRcv_sep' num2str(pairSep) '.mat'];
    save(filename, 'PCs', 'potInterpCosts', 'linInterpCosts', 'xi', 'pairSep', 'pairInds');
    disp(['<-- Saved ' filename]);
end
toc;

%% MOVING RCV, STATIC SRC

cvx_clear;
clear PCs tempPC % just clear a few structs
close all;

disp('******** MOVING RCV, STATIC SRC ********');
pathout = '/Users/aarongeldert/Documents/MATLAB/THESIS/Exported Results/IS_movingRcv_staticSrc/';

roomDims = [7.85 5.35 3.15]; % shoebox room dims from ARTSRAM
roomCen = 0.5*roomDims;
roomVertices = [0,              0; 
                roomDims(1),    0;
                roomDims(1),    roomDims(2);
                0,              roomDims(2)];

srcPos = [6 3.6 1.6]; % arbitrary pos
interpDist = 0.1; % interval between positions
numRcvPos = 26;

rcvPos = ones(numRcvPos, 3) .* [2 4 1.8]; % starting position
rcvDir = [1.5 -1 0]; % direction of rcv movement
rcvDir = interpDist .* rcvDir./vecnorm(rcvDir,2);
rcvPos = rcvPos + (0:(numRcvPos-1)).' * rcvDir;

maxOrder = 10;
PCs = struct('n',cell(1,numRcvPos), 'pos',cell(1,numRcvPos), 'mass',cell(1,numRcvPos));

% error parameters
posError = 0.2;
magError = 0.2;

% downselection target
nTarget = 64;
nRange = 3;

for ii = 1:numRcvPos
    tempPC = imageSourcePointCloud(roomDims, srcPos-roomCen, maxOrder);
    
    % add error
    tempPC.mass = tempPC.mass .* ((1-magError/2)+ magError*rand(size(tempPC.mass)));
    tempPC.pos =  tempPC.pos .* ((1-posError/2)+ posError*rand(size(tempPC.pos)));

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

% interpolation between pairs
for pairSep = 2:25
    toc;
    disp(['-> PAIR SEPARATION: ' num2str(pairSep)]);
    ptsPerPair = pairSep - 1;
    numPairs = numRcvPos - pairSep;
    numPts = ptsPerPair * numPairs;
    pairInds = [1:numPairs; (pairSep+1):numRcvPos].';
    potInterpCosts = NaN(numPairs, ptsPerPair);
    linInterpCosts = NaN(numPairs, ptsPerPair); 

    pairDist = pairSep * interpDist;
    xi = mapTol*pairDist; % OT parameter

    scale = 100;
    offset = 10;
    figure(1); 
    for ii = 1:numPairs 
        pc1 = PCs(pairInds(ii, 1));
        pc2 = PCs(pairInds(ii, 2));
        disp(['Mapping pair ' num2str(ii)]);
        OT = PCOT(pc1, pc2, xi);

        % plot the 2 PCs, room geometry, source and reciever positions 
        cla;
        relSrcPos = srcPos - rcvPos(ii,:);

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

        parfor jj = 1:ptsPerPair

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
            potEval = PCOT(PCot, pcTrue, xiEval);
            linEval = PCOT(PClin, pcTrue, xiEval);

            potInterpCosts(ii,jj) = sum(potEval.Cx.*potEval.Tx,'all');
            linInterpCosts(ii,jj) = sum(linEval.Cx.*linEval.Tx,'all');
        end
    end

    filename = [pathout 'IS_movingRcv_staticSrc_sep' num2str(pairSep) '.mat'];
    save(filename, 'PCs', 'potInterpCosts', 'linInterpCosts', 'xi', 'pairSep', 'pairInds');
    disp(['<-- Saved ' filename]);
end
toc;
%% Combining two PCs
% scales the mass down as well
function PCall = combinePCs(pc1, pc2, k)
    PCall.n = pc1.n + pc2.n;
    PCall.pos = [pc1.pos; pc2.pos];
    PCall.mass = [(1-k)*pc1.mass; k*pc2.mass];
end



