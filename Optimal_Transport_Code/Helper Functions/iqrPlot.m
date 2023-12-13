function [hMedian, hFill] = iqrPlot(x, y, type)
% Created by: Aaron Geldert
% Last modified: 17 Oct 2022

    if nargin < 3
        type = 1;
    end
%     if nargin < 4
%         indB = 0;
%     end
    
%     % I think this doesn't work; we should instead recalculate E 
%     if indB ~= 0
%         y = 10*log10(1+y);
%     end
    
    % y contains distribution data of size [length(x), nPts]
    yLowerQ = prctile(y, 25, 2);
    yMedian = prctile(y, 50, 2);
    yUpperQ = prctile(y, 75, 2);
    
%     yMean = mean(y,2);
%     yLowStd = yMean - std(y,0,2);
%     yUppStd = yMean + std(y,0,2);
    
    xI = linspace(min(x), max(x), 200);

    yLI = pchip(x, yLowerQ, xI);
    yMI = pchip(x, yMedian, xI);
    yUI = pchip(x, yUpperQ, xI);
    
%     yLI = pchip(x, yLowStd, xI);
%     yMI = pchip(x, yMean, xI);
%     yUI = pchip(x, yUppStd, xI);
    
    xPlot = [xI(1:end-1), fliplr(xI)];
    yPlot = [yLI(1:end-1), fliplr(yUI)];
    
    switch type
        case 'lin'
            hFill = fill(xPlot, yPlot, [0.6350 0.0780 0.1840], 'FaceAlpha', 0.15, 'EdgeAlpha',0);
            hold on;
            hMedian = plot(xI, yMI, 'Color', [0.6350 0.0780 0.1840], 'LineStyle', '-');
            hMedian.LineWidth = 1.5;
            hatchfill2(hFill, 'single','HatchAngle',60,'HatchColor',[.7 .1 .1]);
        case 'al'
            hFill = fill(xPlot, yPlot, [0.8500 0.3250 0.0980], 'FaceAlpha', 0.15, 'EdgeAlpha',0);
            hold on;
            hMedian = plot(xI, yMI, 'Color', [0.8500 0.3250 0.0980], 'LineStyle', '-');
            hMedian.LineWidth = 1.5;
            hatchfill2(hFill, 'single','HatchAngle',60,'HatchColor',[.7 .3 .1]);
        case 'nn'
            hFill = fill(xPlot, yPlot, [0.4660 0.6740 0.1880], 'FaceAlpha', 0.15, 'EdgeAlpha',0);
            hold on;
            hMedian = plot(xI, yMI, 'Color', [0.4660 0.6740 0.1880], 'LineStyle', '-');
            hMedian.LineWidth = 2;
            hatchfill2(hFill, 'single','HatchAngle',60,'HatchColor',[.1 .7 .1]);

        case 'ot'
            hFill = fill(xPlot, yPlot, [0 0.4470 0.7410], 'FaceAlpha', 0.15, 'EdgeAlpha',0);
            hold on;
            hMedian = plot(xI, yMI, 'Color', [0 0.4470 0.7410], 'LineStyle', '-');
            hMedian.LineWidth = 2.5;
            hatchfill2(hFill, 'single','HatchAngle',120,'HatchColor',[.1 .1 .7]);
        otherwise
            warning('Unrecognized argument for type.');
    end
    
    
end