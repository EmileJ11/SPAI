function xLs = encodePCtoLS(pc, lsPos, fs)
%ENCODEPCTOLS Takes in a PC and produces loudspeaker signals
%   pc:     point cloud struct with members n, mass, pos
%   lsPos:  loudspeaker positions [L x 3]
%   fs:     sample rate of loudspeaker signals
% Created by: Aaron Geldert
% Last modified: 7 June 2022

c = 343;

% get sample values from distance
sampInds = vecnorm(pc.pos, 2, 2) * fs / c; 
sampInds = round(sampInds) + 1;

% make empty audio signals
xLen = max(sampInds+1); 
numLs = size(lsPos, 1); 
xLs = zeros(xLen, numLs);

% iterate for each point in pc
for ii = 1:pc.n
    
    % find nearest loudspeaker by angle
    angleDist = acosd(pc.pos(ii, :)*lsPos'); 
    [~,lsNear] = min(angleDist);
    
    % assign pressure
    xLs(sampInds(ii),lsNear) = xLs(sampInds(ii),lsNear) + pc.mass(ii);
end

end