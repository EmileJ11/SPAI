function [pos, p] = fillPC(pc, fs, len)
%FILLPC Takes in a PC and produces full vector of DOA and IRs
% DEPRECATED - use encodePCtoLS instead! 7 Jun 2022
% Created by: Aaron Geldert
% Last modified: 31 May 2022

c = 343;

% initialize
pos = zeros(len, 3); % no position if no point
p = zeros(len, 1); % no mass if no point

[sortMass, sortInd] = sort(pc.mass, 'ascend');
sortPos = pc.pos(sortInd,:);

% use distances of each point to get nearest sample index
distInds = vecnorm(sortPos, 2, 2) * fs / c;
% distInds = round(distInds) + 1;

% just in case...
distInds(distInds>len) = len;
distInds(distInds<1) = 1;

for ii = 1:numel(sortMass)
    
    % split up pressure from energy mass
    r = mod(distInds(ii), 1);
    
    % Q: where should the sqrt live? needs to split the energy!
    p(floor(distInds(ii))) = p(floor(distInds(ii))) + (1-r)*sortMass(ii);
    p(ceil(distInds(ii))) = p(ceil(distInds(ii))) + r*sortMass(ii);

    % put the data in the right indices!
    pos(floor(distInds(ii)),:) = sortPos(ii,:); % largest mass determines pos
    pos(ceil(distInds(ii)),:) = sortPos(ii,:);
end

end

