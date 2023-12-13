function SDM = processSDM(srir, fs, filt, threshold, num)
%PROCESSSDM Converts B-format SRIR to full SDM point cloud
% Created by: Aaron Geldert and Nils Meyer-Kahlen
% Last modified: 19 May 2022

if nargin < 3
    filt = 1;
end
if nargin < 4
    threshold = 22;
end
if nargin < 5
    num = 100; % TODO implement a limit on num points
end

c = 343;
n = length(srir);
sampVec = (0:(n-1)).';
tVec = sampVec./fs;
distVec = c*tVec;

omni = srir(:,1); 

% LPF
if filt ~= 0
    [B, A] = butter(4, 5200/(fs/2));
    srirFilt = filtfilt(B, A, srir);
else
    srirFilt = srir;
end

pos = srirFilt(:, 1) .* [srirFilt(:, 4), srirFilt(:, 2), srirFilt(:, 3)];
pos = pos./vecnorm(pos,2,2);
pos = pos.*distVec;

% db thresholding
dbMass = db(omni);
threshInds = find(dbMass > -threshold); 
if(isempty(threshInds))
    warning('No points found above threshold!');
    SDMout = [];
    return;
end

% SDM is a struct with 5 members
SDM.fs = fs; % scalar
SDM.tVec = tVec; % full vector of all possible points
SDM.ir = omni;

SDM.pos = pos(threshInds,:); % cartesian coordinates
SDM.dist = distVec(threshInds); % per point
SDM.mass = dbMass(threshInds)+threshold; % per point
SDM.inds = threshInds;
SDM.n = numel(threshInds);

end

