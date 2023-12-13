function [enVec, tvec] = rirPlot(PC, fs, maxT)
% Created by: Aaron Geldert
% Last modified: 29 Sept 2022
    
    % Takes a point cloud struct and makes mono RIR stem plot
        
    dist = vecnorm(PC.pos,2,2);
    distSamp = round(dist * fs / 343);
    
    numSamp = round(maxT * fs);
    ir = zeros(1,10*numSamp);
    
    for ii = 1:PC.n
        ir(distSamp(ii)) = ir(distSamp(ii)) + PC.mass(ii);
    end
    
    disp(length(ir));
    intTime = 0.004; % 2 ms integration time on each side
    convLen = round(intTime*fs)+1; % odd number for peak symmetry
    
    wind = spreadingFunction(convLen, -20);
    % spreadingFunction returns the same as:
%     wind = 22*(triang(convLen)-1); % triangular masking curve, in dB
%     wind = 10.^(wind./20); % convert to lin
%     wind = wind .* hann(convLen); % smooth the edges
    ir = conv(ir, wind);
    ir = ir(ceil(convLen/2):end);
    disp(length(ir));
    
    sampvec = 1:length(ir);
    tvec = (sampvec-1)./fs;
    plot(tvec, ir, 'Marker', '.');
    xlim([0 maxT]);
    ylim([-0.1 0.6]);
    grid on;
    
    enVec = ir;
    
end