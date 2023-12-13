function SDMout = decimateSDM(SDMin, threshold, num)
%DECIMATESDM Applies a dB threshold to reduce number of points
% Created by: Aaron Geldert and Nils Meyer-Kahlen
% Last modified: 24 May 2022

%% OUTDATED - moved into processSDM.m

% TODO: preprocess via clustering? Downsampling?

if nargin < 3
    num = Inf;
end

% db threshold method
dbMass = db(SDMin.mass);
threshInds = find(dbMass > -threshold); 
if(isempty(threshInds))
    warning('No points found above threshold!');
    SDMout = [];
    return;
end

% TODO: limit it to maximum of num points



SDMout.fs = SDMin.fs;
SDMout.tVec = SDMin.tVec; % full vector
SDMout.pos = SDMin.pos(threshInds, :);
SDMout.dist = SDMin.dist(threshInds, :);
SDMout.mass = SDMin.mass(threshInds, :);

% % convert mass to dB above threshold!
% SDMout.mass = db(SDMin.mass(threshInds, :)) + threshold;

SDMout.inds = threshInds;
SDMout.n = numel(threshInds);

end

