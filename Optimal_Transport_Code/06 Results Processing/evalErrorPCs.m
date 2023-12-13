%% Evaluate the error from different interpolation methods
% Created by: Aaron Geldert
% Last modified: 17 Oct 2022

clear; close all;
folder = '/Users/aarongeldert/Documents/MATLAB/THESIS/ISM_Data/';

errpathout = '/Users/aarongeldert/Documents/MATLAB/THESIS/ISM_Errors/';
% figpathout = '/Volumes/AARON SSD 3/ISM_Figures_Clean/';
figpathout = '/Users/aarongeldert/Documents/MATLAB/THESIS/Fiugres/';


% initialize the error results
LD = load('/Users/aarongeldert/Documents/MATLAB/THESIS/ISM_Data/PC_trapezoidal_order3_dist1_path0.mat');
numPaths = 50;
numk = length(LD.kVec);
ERR.otc = NaN(numk, numPaths);
ERR.otn = NaN(numk, numPaths);
ERR.lic = NaN(numk, numPaths);
ERR.lin = NaN(numk, numPaths);

addpath(folder);
files = dir(folder);
for fileInd = 1:numel(files)
    [pathin,namein,ext] = fileparts(files(fileInd).name);
    if strcmp(ext,'.mat')
        % Each .mat file corresponds to one path
        filename=files(fileInd).name;
        LD = load(filename);
        disp(['Loaded ' filename]);
        
        pathInd = LD.nPath+1;
        kVec = LD.kVec;
        for kInd = 1:(length(kVec))
            k = kVec(kInd);
            if k~=0 && k~=1
                % intermediate point
                ERR.otc(kInd, pathInd) = getIRerror(LD.PCotc(kInd), LD.PCgtc(kInd), LD.fs);
                ERR.otn(kInd, pathInd) = getIRerror(LD.PCotn(kInd), LD.PCgtc(kInd), LD.fs); % clean ground truth always!
                ERR.lic(kInd, pathInd) = getIRerror(LD.PClic(kInd), LD.PCgtc(kInd), LD.fs);
                ERR.lin(kInd, pathInd) = getIRerror(LD.PClin(kInd), LD.PCgtc(kInd), LD.fs);
            else
                % endpoint
                ERR.otc(kInd, pathInd) = 0;
                ERR.otn(kInd, pathInd) = 0;
                ERR.lic(kInd, pathInd) = 0;
                ERR.lin(kInd, pathInd) = 0;
            end
        end
    end
end

%%
close all;
figure(1);
% plot(kVec, median(ERR.lic,2));
% hold on;
% plot(kVec, median(ERR.lin,2));
% plot(kVec, median(ERR.otc,2));
% plot(kVec, median(ERR.otn,2));
% legend('Lin Clean', 'Lin Noisy', 'POT Clean', 'POT Noisy');

[hMlic, hFlic] = iqrPlot(kVec, ERR.lic);
[hMotc, hFotc] = iqrPlot(kVec, ERR.otc);

