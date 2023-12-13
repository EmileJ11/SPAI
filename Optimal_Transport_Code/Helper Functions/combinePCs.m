function PCall = combinePCs(pc1, pc2, k)
%COMBINEPCS Implements a simple linear combination of the point clouds
%where k is the interpolation parameter.
% e.g. PCall is pc1 when k=0, and pc2 when k=1
% Resulting PC may contain multiple masses in the same positions.
% Created by: Aaron Geldert
% Last modified: 5 Oct 2022
    PCall.n = pc1.n + pc2.n;
    PCall.pos = [pc1.pos; pc2.pos];
    PCall.mass = [(1-k)*pc1.mass; k*pc2.mass];
end