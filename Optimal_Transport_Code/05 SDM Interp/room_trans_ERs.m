%% Turning the SDM / Diffuse datasets into clustered PCs
% Created by: Aaron Geldert
% Last modified: 21 Aug 2022

clear; clc; close all;

folders = ["/Users/aarongeldert/Documents/MATLAB/THESIS/Exported Results/RoomTrans_EnergyDiffuse/"];

% Output .mat file location
pathout = '/Users/aarongeldert/Documents/MATLAB/THESIS/Exported Results/RoomTrans_PCs/';

c = 343;

% Early reflection parameters
p_ert = 0.06; % arrival time threshold
p_dbt = -60; % energy x prominence threshold (dB)

% DOA averaging window (matches analysis window)
winLen = 32; % window length
halfLen = round(winLen/2);
win = hamming(winLen); % window

% Load the data
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
            M = size(LD.energy,2);
            
            assert(~any(isnan(LD.doa),'all'), ['NaN found in DOA data in ' namein]);
            disp('Loaded');
            
            % v1) multiply the diffuse and energy vectors (?)
            erVec = LD.diffuse.*LD.energy;
            
            % time window
            maxSamp = p_ert*LD.fs;
            erVec = erVec(1:maxSamp,:);
            
            % preallocate PCs
            allPCs = struct('n', cell(1, M), 'pos', cell(1, M), 'mass', cell(1, M));
            
            for measInd = 1:M
                
                % maybe renormalize for each position...?
                erVec(:,measInd) = erVec(:,measInd) ./ max(erVec(:,measInd));
                
                [pks, pkInds] = findpeaks(db(erVec(:,measInd)), 'MinPeakHeight', p_dbt);
                % [sortedPks, sortInds] = sort(pks,'descend');
                % maybe have a sorting of largest peaks perhaps...?
                
                numPeaks = numel(pkInds);
                
                % DOA data per peak
                pos = NaN(numPeaks, 3); % cartesian data only!
                
                % iterate through each peak and find window avg. DOA
                for pki = 1:numPeaks
                    
                    % v1: weighted average DOA in window
%                     winSamps = (pkInds(pki)-halfLen):(pkInds(pki)+halfLen-1);
%                     valSamps = winSamps>0 & winSamps<=maxSamp; % logical
%                     posVals = zeros(winLen,3);
%                     posVals(valSamps,:) = S2C(LD.doa(winSamps(valSamps),:,measInd));
%                     posVals = posVals .* win;
%                     dir = mean(posVals,1);
                    
                    % v2: just at the peak of directionality * energy
                    dir = S2C(LD.doa(pkInds(pki), :, measInd));
                    
                    dir = dir./vecnorm(dir,2,2); % normalize to unit vector
                    dist = LD.tvec(pkInds(pki)) * c + LD.dist(measInd); % delay time and src-rcv correction 
                    pos(pki,:) = dir .* dist;
                     
                end
                
                % make a point cloud data struct
                allPCs(measInd).pos = pos;
                allPCs(measInd).mass = sqrt(LD.energy(pkInds,measInd)); % sqrt converts energy to pressure
                allPCs(measInd).n = length(pkInds);
                disp(['M: ' num2str(measInd) ' num pts: ' num2str(allPCs(measInd).n) ]);
            end
            
            % save .mat data
            fileOutName = [pathout 'PCs_' namein '.mat'];
            save(fileOutName, 'allPCs', 'p_ert', 'p_dbt');
            disp(['Saved ' fileOutName]);
            
        end
    end
end
toc;

%% Visualize the PCs of the measurements

clear; clc; close all;

folders = ["/Users/aarongeldert/Documents/MATLAB/THESIS/Exported Results/RoomTrans_PCs/"];

% Output .gif file location
pathout = '/Users/aarongeldert/Documents/MATLAB/THESIS/Exported Results/SDM Room Transitions/';

figure(1);
axlim = [-16 16 -16 16 -16 16];
scale = 80;
c = 343;
fs = 48e3;

tgif = 16; % gif time, seconds

tdsgn = getTdesign(6);

% Load the data
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
            
            % prepare to generate a GIF
            gifname_SDM = [pathout 'SDM_ER_' cleanName, '.gif'];
            gifname_LS = [pathout 'LS_ER_' cleanName, '.gif'];
            
            for ii = 1:M
                
                figure(1); cla;
                scatter3(0,0,0,'*'); % origin  (listener position)
                hold on; grid on;
                scatter3(PC(ii).pos(:,1), PC(ii).pos(:,2), PC(ii).pos(:,3),...
                    PC(ii).mass * scale, 'ko','filled');
                view([90 90]);
                xlabel('X'); ylabel('Y'); zlabel('Z');
                axis(axlim);
                title([strrep(cleanName, '_',' ') newline 'pos ' num2str(ii) ' / ' num2str(M)]);
                
                % save point cloud pic to the gif
                frame = getframe(gcf); 
                im = frame2im(frame); 
                [imind,cm] = rgb2ind(im,256); 
                if ii==1 
                    imwrite(imind,cm,gifname_SDM,'gif', 'Loopcount',inf,'DelayTime',tgif/M); 
                else 
                    imwrite(imind,cm,gifname_SDM,'gif','WriteMode','append','DelayTime',tgif/M); 
                end
                
                figure(2); cla;
                LS = 2 * encodePCtoLS(PC(ii), tdsgn, fs);
                plot(LS + (0:23)); xlim([1 LD.p_ert*fs]); 
                ylim([0 25]);
                
                % save ls pic to the gif
                frame = getframe(gcf); 
                im = frame2im(frame); 
                [imind,cm] = rgb2ind(im,256); 
                if ii==1 
                    imwrite(imind,cm,gifname_LS,'gif', 'Loopcount',inf,'DelayTime',tgif/M); 
                else 
                    imwrite(imind,cm,gifname_LS,'gif','WriteMode','append','DelayTime',tgif/M); 
                end
                
            end
            toc;
        end
    end
end

toc;
