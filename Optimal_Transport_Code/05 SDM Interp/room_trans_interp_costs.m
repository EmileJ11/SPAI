%% Cost of Optimal Transport of PCs from Room Trans dataset
% Created by: Aaron Geldert
% Last modified: 21 Aug 2022

% 1) Load the PCs from the dataset
clear; clc; close all;
folders = ["/Users/aarongeldert/Documents/MATLAB/THESIS/Exported Results/RoomTrans_PCs/"];

interpDist = 0.05; % 5 cm dist specified
xi = 4*interpDist; % OT parameter

% Output .mat file location
pathout = '/Users/aarongeldert/Documents/MATLAB/THESIS/Exported Results/RoomTrans_InterpCosts/';

tic;
for ff = 1:numel(folders)
    addpath(folders(ff));
    files=dir(folders(ff));
    
    for nn = 1:numel(files)
        [~,namein,ext] = fileparts(files(nn).name);
        if isempty(namein)
            continue;
        elseif strcmp(namein(1),'.')
            continue;
        end

        if strcmp(ext,'.mat')
            % load data
            filename = [namein ext];
            LD = load(filename);
            cleanName = erase(namein, 'PCs_analyses_');
            
            PC = LD.allPCs;
            M = length(PC);
            
            % 2) OPTIMAL TRANSPORT
            % Compute OT mapping between a pair of measurement positions
            
            pairSep = 2; % interp pos 2 from 1 and 3
            ptsPerPair = pairSep - 1;
            numPairs = M - pairSep;
            numPts = ptsPerPair * numPairs;
            pairInds = [1:numPairs; (pairSep+1):M].';
            potInterpCosts = NaN(numPairs, ptsPerPair);
            linInterpCosts = NaN(numPairs, ptsPerPair);
            
            for ii = 1:numPairs
                
                pc1 = PC(pairInds(ii, 1));
                pc2 = PC(pairInds(ii, 2));
                OT = PCOT(pc1, pc2, xi);
                
                for jj = 1:ptsPerPair
                    
                    % 3) INTERPOLATION
                    % a) Create OT interpolated PC
                    k = jj./pairSep;
                    [PCot, ~, ~] = OT.interpPC(k);
                    
                    % b) Create Linear combination PC
                    PClin = combinePCs(pc1, pc2);
                    
                    % c) Create NN interpolated PC (?)
                    
                    % 4) OBJECTIVE EVALUATION
                    % Compute OT cost between interpolated PC and ground truth PC
                    % use no virtual masses (transport entire point clouds)
                    
                    xiEval = xi; % same number (why? dunno!)
                    pcTrue = PC(pairInds(ii, 1) + jj);
                    
                    potEval = PCOT(PCot, pcTrue, xiEval);
                    linEval = PCOT(PClin, pcTrue, xiEval);
                    
                    potInterpCosts(ii,jj) = sum(potEval.Cx.*potEval.Tx,'all');
                    linInterpCosts(ii,jj) = sum(linEval.Cx.*linEval.Tx,'all');
                    
                end
                
                disp(['Computing pair ' num2str(ii)]);
                
            end

            fileOutName = [pathout 'interp_costs_' cleanName '.mat'];
            save(fileOutName, 'potInterpCosts', 'linInterpCosts', 'xi', 'pairSep', 'pairInds');
            
        end
    end
end
toc;

function PClin = combinePCs(pc1, pc2)
    % Linear Combination of PCs
    PClin.pos = [pc1.pos; pc2.pos];
    PClin.mass = 0.5 * [pc1.mass; pc2.mass];
    PClin.n = pc1.n + pc2.n;
end
