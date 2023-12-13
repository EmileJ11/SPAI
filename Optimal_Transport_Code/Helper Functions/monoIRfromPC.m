function [ir, tvec] = monoIRfromPC(PC, maxT, fs)
% maxT is IR length in seconds
% Created by: Aaron Geldert
% Last modified: 10 Oct 2022
    
    if nargin<3
        fs = 48e3;
    end
    
    % Takes a point cloud struct and makes mono RIR
    dist = vecnorm(PC.pos,2,2);
    distSamp = round(dist * fs / 343);

    numSamp = round(maxT * fs);
    
    ir = zeros(1,max(max(distSamp), numSamp));
    for ii = 1:PC.n
        ir(distSamp(ii)) = ir(distSamp(ii)) + PC.mass(ii);
    end
    
    ir = ir(1:numSamp);
    tvec = linspace(0,maxT,length(ir));
end