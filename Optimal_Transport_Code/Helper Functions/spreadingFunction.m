function sf = spreadingFunction(convLen, dbExtent, pl)
%SPREADINGFUNCTION produces a window that approximates a temporal masking 
% curve.
% Created by: Aaron Geldert
% Last modified: 4 Oct 2022

if nargin < 3
    pl = 0;
end
if nargin < 2
    % defines slope for -20 dB (x0.1) at window edges
    dbExtent = -10; % negative for clarity, later multiplied by -1 again
end

convLen = round(convLen);
halfLen = floor(convLen/2);
qrtrLen = floor(halfLen/2);

sf = -dbExtent*(triang(convLen)-1); % triangular masking curve, in dB
sf = 10.^(sf./20); % convert to lin

% smooth the window edges to zero
% halfHann = hann(halfLen);
% smoothWin = [halfHann(1:qrtrLen); ones(halfLen+1,1); halfHann(qrtrLen+1:halfLen)];
% sf = sf .* smoothWin;
sf = sf .* sqrt(hann(convLen));

if pl ~= 0
    % Plots for the thesis:
    fs = 48e3;
    halfLen = floor(convLen/2);
    tvec = (-halfLen:halfLen)./(fs/1e3); % in ms
    figure(2);
    subplot(121);
    plot(tvec, sf); % linear
    axis([min(tvec) max(tvec) -0.05 1.05]);
    ylabel('Amplitude');
    xlabel('Time (ms)');
    xticks(-2:2)
    grid on;
    yline(0); xline(0);
    set(gca, 'fontName', 'Times')

    subplot(122)
    plot(tvec, db(sf)); % in dB
    axis([min(tvec) max(tvec) -25 1]);
    ylabel('dB');
    xlabel('Time (ms)');
    xticks(-2:2)
    grid on;
    yline(0); xline(0);
    set(gca, 'fontName', 'Times')
    sgtitle('Spreading function', 'fontName', 'Times')
end
end