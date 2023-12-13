%% Apply diffuseness mask to energy
% Created by: Aaron Geldert
% Last modified: 23 Aug 2022

clear; clc; close all;

folders = ["/Users/aarongeldert/Documents/MATLAB/THESIS/Exported Results/RoomTrans_EnergyDiffuse/"];

% Output .gif file location
pathout = '/Users/aarongeldert/Documents/MATLAB/THESIS/Exported Results/DiffuseEstGifs/';

tgif = 16; % gif time, seconds
c = 343;

% Early reflection parameters
p_ert = 0.06; % arrival time threshold
d_thres = 0.16; % diffuse threshold

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
            disp('Loaded');
            
            % get diffuseness mask
            isDirectional = LD.diffuse > d_thres;
            
            erVec = LD.energy;
            erVec = erVec ./ max(erVec, [], 'all'); % normalize
            erVec = erVec .* isDirectional; % apply mask
            
            % time window
            maxSamp = round(p_ert*LD.fs);
            
            % prepare to generate a GIF
            cleanName = erase(namein, 'analyses_');
            gifname = [pathout cleanName '.gif'];

            for measInd = 1:M
                
                figure(1); cla;
                
                plot(erVec(:, measInd));
                axis([0 maxSamp 0 1]);
                title([strrep(cleanName,'_',' ') newline 'pos ' num2str(measInd) ' / ' num2str(M)])
                xlabel('Time (samples)');
                ylabel('Energy (p^2)');
                
                % save ls pic to the gif
                frame = getframe(gcf); 
                im = frame2im(frame); 
                [imind,cm] = rgb2ind(im,256); 
                if measInd==1 
                    imwrite(imind,cm,gifname,'gif', 'Loopcount',inf,'DelayTime',tgif/M); 
                else 
                    imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',tgif/M); 
                end
            end
                
        end
    end
end
toc;

