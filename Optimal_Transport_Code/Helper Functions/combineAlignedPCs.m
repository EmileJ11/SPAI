function PCall = combineAlignedPCs(pc1, pc2, k)
%COMBINEALIGNEDPCS Implements a simple linear combination of the point clouds
%after direct sound alignment, where k is the interpolation parameter.
% e.g. PCall is pc1 when k=0, and pc2 when k=1
% Resulting PC may contain multiple masses in the same positions.
% Created by: Aaron Geldert
% Last modified: 23 Nov 2022

    dist1 = vecnorm(pc1.pos,2,2); % distance (m) of each point
    dist2 = vecnorm(pc2.pos,2,2);

    directDist1 = min(dist1); % distance asspciated with direct sound (TOA)
    directDist2 = min(dist2);
    directDistk = (1-k)*directDist1 + k*directDist2; % expected interpolated direct dist
    
    alignedDist1 = (dist1-directDist1)+directDistk;
    alignedDist2 = (dist2-directDist2)+directDistk;
    
    alignedPos1 = alignedDist1.*(pc1.pos./dist1);
    alignedPos2 = alignedDist2.*(pc2.pos./dist2);

    PCall.n = pc1.n + pc2.n;
    PCall.pos = [alignedPos1; alignedPos2];
    PCall.mass = [(1-k)*pc1.mass; k*pc2.mass];
end