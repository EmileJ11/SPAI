function [doas, p] = fillPointCloud(SDM)
%FILLPOINTCLOUD Takes in a SDM point cloud and produces full vector of IRs
% Created by: Aaron Geldert
% Last modified: 19 May 2022

len = length(SDM.tVec); % full time vector
c = 343;

% initialize
doas = NaN(len, 3); % no position if no point
p = zeros(len, 1); % no mass if no point

% use distances of each point to get nearest sample index
distInds = SDM.dist * SDM.fs / c;
distInds = round(distInds) + 1;

% just in case...
distInds(distInds>len) = len;
distInds(distInds<1) = 1;

% put the data in the right indices!
doas(distInds,:) = SDM.pos;
p(distInds) = SDM.mass;

end

