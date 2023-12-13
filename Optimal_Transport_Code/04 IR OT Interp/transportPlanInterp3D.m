function [DOA, M] = transportPlanInterp3D(k, doa1, doa2, T)
%TRANSPORTPLANINTERP Summary of this function goes here
%   Interpolate at k, using transport plans T and T2
% Parameters:
%   k: interpolation parameter (0:1)
%   doa1:
% 
%   T: transport plan (masses placed places, zero otherwise)
% Created by: Aaron Geldert
% Last Modified: 2 May 2022

    assert(isequal(size(doa1),size(doa2)), 'Doa sizes must match.');
    assert(k>=0 && k<=1, 'k must be between 0 and 1.');

    % distLen = length(T); % this length does not matter anymore
    [I, J, M] = find(T);
    numMasses = numel(M);
    assert(length(doa1)==numMasses,'DOA vector lengths must match the number of nonzero elements in T.')

    % get the interpolated dist from X->Y
    DOA = zeros(numMasses,3);

    % mass locations are no longer quantized into bins!
    for ii = 1:numMasses

        DOA(ii,:) = (1-k).*doa1(I(ii),:) + k.*doa2(J(ii),:);

    end

end

